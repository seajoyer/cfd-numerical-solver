#include "output/GIFWriter.hpp"

#include <vtkAxis.h>
#include <vtkChartLegend.h>
#include <vtkChartXY.h>
#include <vtkContextScene.h>
#include <vtkContextView.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkPen.h>
#include <vtkPlot.h>
#include <vtkPlotLine.h>
#include <vtkPointData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkTextProperty.h>
#include <vtkUnsignedCharArray.h>
#include <vtkWindowToImageFilter.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "data/DataLayer.hpp"
#include "utils/StringUtils.hpp"

// Include gif.h for GIF encoding
// This is a single-header library: https://github.com/charlietangora/gif-h
#include "external/gif.h"

/**
 * @brief Frame data structure holding raw RGBA pixels
 */
struct Frame {
    std::vector<uint8_t> pixels;  // RGBA data (width * height * 4 bytes)
    int width;
    int height;
};

// PIMPL implementation
class GIFWriter::Impl {
   public:
    mutable std::vector<Frame> frames;  // Mutable to allow modification in const Write()
    
    Impl() = default;
    
    void AddFrame(const std::vector<uint8_t>& pixels, int width, int height) {
        Frame frame;
        frame.pixels = pixels;
        frame.width = width;
        frame.height = height;
        frames.push_back(std::move(frame));
    }
    
    void Clear() {
        frames.clear();
    }
    
    [[nodiscard]] auto GetFrameCount() const -> std::size_t {
        return frames.size();
    }
};

GIFWriter::GIFWriter(const std::string& output_dir, int width, int height,
                     int delay_centiseconds)
    : output_dir_(output_dir),
      width_(width),
      height_(height),
      delay_centiseconds_(delay_centiseconds),
      pimpl_(std::make_unique<Impl>()) {
    // Create output directory if it doesn't exist
    std::filesystem::create_directories(output_dir_);
}

GIFWriter::~GIFWriter() = default;

GIFWriter::GIFWriter(GIFWriter&&) noexcept = default;
auto GIFWriter::operator=(GIFWriter&&) noexcept -> GIFWriter& = default;

auto GIFWriter::GenerateFilename(const Settings& settings) const -> std::string {
    std::ostringstream oss;
    oss << output_dir_ << "/" << settings.solver << "__R_" << settings.reconstruction
        << "__N_" << settings.N << "__CFL_" << utils::DoubleWithoutDot(settings.cfl)
        << ".gif";
    return oss.str();
}

auto GIFWriter::GetFrameCount() const -> std::size_t {
    return pimpl_->GetFrameCount();
}

void GIFWriter::SetFrameDelay(int centiseconds) {
    delay_centiseconds_ = centiseconds;
}

auto GIFWriter::GetFrameDelay() const -> int {
    return delay_centiseconds_;
}

void GIFWriter::Write(const DataLayer& layer, const Settings& settings,
                      std::size_t step, double time) const {
    if (layer.GetDim() >= 2) {
        static bool warned = false;
        if (!warned) {
            std::cerr << "GIFWriter: 2D visualization not yet supported, skipping GIF output.\n"
                      << "  Use VTK format for 2D simulations.\n";
            warned = true;
        }
        return;
    }

    Write(layer, nullptr, settings, step, time);
}

void GIFWriter::Write(const DataLayer& layer, const DataLayer* analytical_layer,
                      const Settings& settings, std::size_t step, double time) const {
    if (layer.GetDim() >= 2) {
        return;  // Warning already issued
    }
    RenderFrame(layer, analytical_layer, settings, step, time);
}

void GIFWriter::RenderFrame(const DataLayer& layer, const DataLayer* analytical_layer,
                            const Settings& settings, std::size_t step,
                            double time) const {
    const int start = layer.GetCoreStart();
    const int end = layer.GetCoreEndExclusive();
    const int n_points = end - start;

    if (n_points <= 0) {
        throw std::runtime_error("Invalid core range for GIF frame");
    }

    // Calculate scale factor relative to baseline 1200x900
    const float baseline_width = 1200.0f;
    const float baseline_height = 900.0f;
    const float width_ratio = static_cast<float>(width_) / baseline_width;
    const float height_ratio = static_cast<float>(height_) / baseline_height;
    const float scale_factor = std::sqrt(width_ratio * height_ratio);

    // Baseline values optimized for 1200x900
    const float base_title_font_size = 18.0f;
    const float base_axis_title_font_size = 18.0f;
    const float base_label_font_size = 14.0f;
    const float base_line_width = 2.5f;
    const float base_margin = 10.0f;

    // Apply scaling
    const int title_font_size =
        static_cast<int>(base_title_font_size * scale_factor + 0.5f);
    const int axis_title_font_size =
        static_cast<int>(base_axis_title_font_size * scale_factor + 0.5f);
    const int label_font_size =
        static_cast<int>(base_label_font_size * scale_factor + 0.5f);
    const float line_width = base_line_width * scale_factor;
    const float margin = base_margin * scale_factor;

    // Create context view for off-screen rendering
    vtkNew<vtkContextView> view;
    vtkRenderWindow* renderWindow = view->GetRenderWindow();
    renderWindow->SetOffScreenRendering(1);
    renderWindow->SetSize(width_, height_);
    view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

    // Create 4 charts
    vtkNew<vtkChartXY> chart_density;
    vtkNew<vtkChartXY> chart_velocity;
    vtkNew<vtkChartXY> chart_pressure;
    vtkNew<vtkChartXY> chart_energy;

    view->GetScene()->AddItem(chart_density);
    view->GetScene()->AddItem(chart_velocity);
    view->GetScene()->AddItem(chart_pressure);
    view->GetScene()->AddItem(chart_energy);

    // Calculate chart dimensions
    const float w = static_cast<float>(width_);
    const float h = static_cast<float>(height_);
    const float chart_w = w / 2.0f;
    const float chart_h = h / 2.0f;

    // Position charts in 2x2 grid
    chart_density->SetAutoSize(false);
    chart_density->SetSize(
        vtkRectf(margin, chart_h + margin, chart_w - 2 * margin, chart_h - 2 * margin));

    chart_velocity->SetAutoSize(false);
    chart_velocity->SetSize(vtkRectf(chart_w + margin, chart_h + margin,
                                     chart_w - 2 * margin, chart_h - 2 * margin));

    chart_pressure->SetAutoSize(false);
    chart_pressure->SetSize(
        vtkRectf(margin, margin, chart_w - 2 * margin, chart_h - 2 * margin));

    chart_energy->SetAutoSize(false);
    chart_energy->SetSize(vtkRectf(chart_w + margin, margin, chart_w - 2 * margin,
                                   chart_h - 2 * margin));

    // Helper lambda to create data table
    auto createTable = [&](const DataLayer& data_layer, const std::string& y_name,
                           auto accessor) -> vtkSmartPointer<vtkTable> {
        vtkNew<vtkTable> table;

        vtkNew<vtkFloatArray> x_array;
        x_array->SetName("x");
        x_array->SetNumberOfValues(n_points);

        vtkNew<vtkFloatArray> y_array;
        y_array->SetName(y_name.c_str());
        y_array->SetNumberOfValues(n_points);

        for (int i = 0; i < n_points; ++i) {
            const int idx = start + i;
            x_array->SetValue(i, static_cast<float>(data_layer.xc(idx)));
            y_array->SetValue(i, static_cast<float>(accessor(data_layer, idx)));
        }

        table->AddColumn(x_array);
        table->AddColumn(y_array);
        return table;
    };

    // Accessors
    auto rho_accessor = [](const DataLayer& dl, int i) { return dl.rho(i); };
    auto u_accessor = [](const DataLayer& dl, int i) { return dl.u(i); };
    auto P_accessor = [](const DataLayer& dl, int i) { return dl.P(i); };
    auto U_accessor = [](const DataLayer& dl, int i) { return dl.U(i); };

    // Create tables
    auto table_rho = createTable(layer, "Density", rho_accessor);
    auto table_u = createTable(layer, "Velocity", u_accessor);
    auto table_P = createTable(layer, "Pressure", P_accessor);
    auto table_U = createTable(layer, "Energy", U_accessor);

    // Configure chart appearance
    auto configureChart = [&](vtkChartXY* chart, const std::string& title,
                              const std::string& y_label) {
        chart->SetTitle(title);
        chart->GetTitleProperties()->SetFontSize(title_font_size);
        chart->GetTitleProperties()->SetBold(true);
        chart->GetTitleProperties()->SetColor(0.0, 0.0, 0.0);

        chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("x");
        chart->GetAxis(vtkAxis::BOTTOM)->GetTitleProperties()->SetFontSize(
            axis_title_font_size);
        chart->GetAxis(vtkAxis::BOTTOM)->GetTitleProperties()->SetColor(0.0, 0.0, 0.0);
        chart->GetAxis(vtkAxis::BOTTOM)->GetLabelProperties()->SetFontSize(
            label_font_size);
        chart->GetAxis(vtkAxis::BOTTOM)->GetLabelProperties()->SetColor(0.0, 0.0, 0.0);
        chart->GetAxis(vtkAxis::BOTTOM)->GetPen()->SetColor(0, 0, 0);
        chart->GetAxis(vtkAxis::BOTTOM)->GetGridPen()->SetColor(200, 200, 200);

        chart->GetAxis(vtkAxis::LEFT)->SetTitle(y_label);
        chart->GetAxis(vtkAxis::LEFT)->GetTitleProperties()->SetFontSize(
            axis_title_font_size);
        chart->GetAxis(vtkAxis::LEFT)->GetTitleProperties()->SetColor(0.0, 0.0, 0.0);
        chart->GetAxis(vtkAxis::LEFT)->GetLabelProperties()->SetFontSize(label_font_size);
        chart->GetAxis(vtkAxis::LEFT)->GetLabelProperties()->SetColor(0.0, 0.0, 0.0);
        chart->GetAxis(vtkAxis::LEFT)->GetPen()->SetColor(0, 0, 0);
        chart->GetAxis(vtkAxis::LEFT)->GetGridPen()->SetColor(200, 200, 200);

        chart->SetShowLegend(true);
        chart->GetLegend()->SetHorizontalAlignment(2);  // RIGHT
        chart->GetLegend()->SetVerticalAlignment(1);    // TOP
        chart->GetLegend()->GetLabelProperties()->SetFontSize(label_font_size);
    };

    // Helper to add plot
    auto addPlot = [&](vtkChartXY* chart, vtkTable* table, const std::string& y_column,
                       unsigned char r, unsigned char g, unsigned char b,
                       const std::string& label) -> void {
        vtkPlot* plot = chart->AddPlot(vtkChart::LINE);
        plot->SetInputData(table, "x", y_column);
        vtkPen* pen = plot->GetPen();
        pen->SetWidth(line_width);
        pen->SetColor(r, g, b, 255);
        plot->SetLabel(label);
    };

    // Configure charts with titles
    std::ostringstream title_oss;
    title_oss << std::fixed << std::setprecision(4);

    title_oss.str("");
    title_oss << "Density (t=" << time << ")";
    configureChart(chart_density, title_oss.str(), "rho");

    title_oss.str("");
    title_oss << "Velocity (t=" << time << ")";
    configureChart(chart_velocity, title_oss.str(), "u");

    title_oss.str("");
    title_oss << "Pressure (t=" << time << ")";
    configureChart(chart_pressure, title_oss.str(), "P");

    title_oss.str("");
    title_oss << "Specific Internal Energy (t=" << time << ")";
    configureChart(chart_energy, title_oss.str(), "e");

    // Add analytical solution first (black line)
    if (analytical_layer != nullptr) {
        auto table_rho_a = createTable(*analytical_layer, "Density", rho_accessor);
        auto table_u_a = createTable(*analytical_layer, "Velocity", u_accessor);
        auto table_P_a = createTable(*analytical_layer, "Pressure", P_accessor);
        auto table_U_a = createTable(*analytical_layer, "Energy", U_accessor);

        addPlot(chart_density, table_rho_a, "Density", 0, 0, 0, "Analytical");
        addPlot(chart_velocity, table_u_a, "Velocity", 0, 0, 0, "Analytical");
        addPlot(chart_pressure, table_P_a, "Pressure", 0, 0, 0, "Analytical");
        addPlot(chart_energy, table_U_a, "Energy", 0, 0, 0, "Analytical");
    }

    // Add numerical solution (red line)
    addPlot(chart_density, table_rho, "Density", 255, 0, 0, "Numerical");
    addPlot(chart_velocity, table_u, "Velocity", 255, 0, 0, "Numerical");
    addPlot(chart_pressure, table_P, "Pressure", 255, 0, 0, "Numerical");
    addPlot(chart_energy, table_U, "Energy", 255, 0, 0, "Numerical");

    // Render
    renderWindow->SetMultiSamples(0);
    renderWindow->SetLineSmoothing(1);
    renderWindow->SetPolygonSmoothing(1);
    view->GetScene()->SetDirty(true);
    renderWindow->Render();
    renderWindow->Render();

    // Capture frame
    vtkNew<vtkWindowToImageFilter> windowToImageFilter;
    windowToImageFilter->SetInput(renderWindow);
    windowToImageFilter->SetScale(1);
    windowToImageFilter->SetInputBufferTypeToRGBA();
    windowToImageFilter->ReadFrontBufferOff();
    windowToImageFilter->Update();

    // Extract pixel data
    vtkImageData* imageData = windowToImageFilter->GetOutput();
    int* dims = imageData->GetDimensions();
    int frame_width = dims[0];
    int frame_height = dims[1];

    vtkUnsignedCharArray* pixels = vtkUnsignedCharArray::SafeDownCast(
        imageData->GetPointData()->GetScalars());

    if (!pixels) {
        throw std::runtime_error("Failed to get pixel data from rendered frame");
    }

    // Copy pixel data to vector (RGBA format)
    std::size_t pixel_count = static_cast<std::size_t>(frame_width * frame_height * 4);
    std::vector<uint8_t> frame_pixels(pixel_count);

    // VTK stores images bottom-up, so we need to flip vertically for GIF
    for (int y = 0; y < frame_height; ++y) {
        int src_row = frame_height - 1 - y;
        for (int x = 0; x < frame_width; ++x) {
            int src_idx = (src_row * frame_width + x) * 4;
            int dst_idx = (y * frame_width + x) * 4;
            frame_pixels[dst_idx + 0] = pixels->GetValue(src_idx + 0);  // R
            frame_pixels[dst_idx + 1] = pixels->GetValue(src_idx + 1);  // G
            frame_pixels[dst_idx + 2] = pixels->GetValue(src_idx + 2);  // B
            frame_pixels[dst_idx + 3] = pixels->GetValue(src_idx + 3);  // A
        }
    }

    pimpl_->AddFrame(frame_pixels, frame_width, frame_height);
}

auto GIFWriter::Finalize(const Settings& settings) -> std::string {
    if (pimpl_->frames.empty()) {
        std::cerr << "Warning: No frames to write to GIF\n";
        return "";
    }

    std::string filename = GenerateFilename(settings);

    // Initialize GIF writer
    GifWriter gif_writer;
    if (!GifBegin(&gif_writer, filename.c_str(), width_, height_,
                  delay_centiseconds_)) {
        throw std::runtime_error("Failed to initialize GIF writer for: " + filename);
    }

    // Write each frame
    for (const auto& frame : pimpl_->frames) {
        // gif.h expects RGBA data with dimensions matching the GIF
        if (frame.width != width_ || frame.height != height_) {
            std::cerr << "Warning: Frame size mismatch, skipping frame\n";
            continue;
        }

        if (!GifWriteFrame(&gif_writer, frame.pixels.data(), frame.width, frame.height,
                           delay_centiseconds_)) {
            std::cerr << "Warning: Failed to write frame to GIF\n";
        }
    }

    // Finalize GIF
    if (!GifEnd(&gif_writer)) {
        throw std::runtime_error("Failed to finalize GIF: " + filename);
    }

    // Clear frame buffer to free memory
    pimpl_->Clear();

    return filename;
}
