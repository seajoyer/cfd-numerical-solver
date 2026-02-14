#include "output/PNGWriter.hpp"

#include <vtkAxis.h>
#include <vtkChartLegend.h>
#include <vtkChartXY.h>
#include <vtkContextScene.h>
#include <vtkContextView.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkPNGWriter.h>
#include <vtkPen.h>
#include <vtkPlot.h>
#include <vtkPlotLine.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkTextProperty.h>
#include <vtkWindowToImageFilter.h>

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "data/DataLayer.hpp"

// PIMPL implementation
class PNGWriter::Impl {
   public:
    // No persistent state needed - we create fresh objects for each write
    Impl() = default;
};

PNGWriter::PNGWriter(const std::string& output_dir, int width, int height)
    : output_dir_(output_dir),
      width_(width),
      height_(height),
      pimpl_(std::make_unique<Impl>()) {
    // Create output directory if it doesn't exist
    std::filesystem::create_directories(output_dir_);
}

PNGWriter::~PNGWriter() = default;

auto PNGWriter::GenerateFilename(std::size_t step) const -> std::string {
    std::ostringstream oss;
    oss << output_dir_ << "/step_" << std::setw(4) << std::setfill('0') << step << ".png";
    return oss.str();
}

void PNGWriter::Write(const DataLayer& layer, const Settings& settings,
                      std::size_t step, double time) const {
    if (layer.GetDim() >= 2) {
        static bool warned = false;
        if (!warned) {
            std::cerr << "PNGWriter: 2D visualization not yet supported, skipping PNG output.\n"
                      << "  Use VTK format for 2D simulations.\n";
            warned = true;
        }
        return;
    }

    Write(layer, nullptr, settings, step, time);
}

void PNGWriter::Write(const DataLayer& layer, const DataLayer* analytical_layer,
                      const Settings& settings, std::size_t step, double time) const {
    if (layer.GetDim() >= 2) {
        return;  // Warning already issued by the other Write overload
    }

    const int start = layer.GetCoreStart();
    const int end = layer.GetCoreEndExclusive();
    const int n_points = end - start;

    if (n_points <= 0) {
        throw std::runtime_error("Invalid core range for PNG output");
    }

    // Calculate scale factor relative to baseline 1200x900
    // Use geometric mean of width and height ratios for balanced scaling
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
    
    // Apply scaling to all visual elements
    const int title_font_size = static_cast<int>(base_title_font_size * scale_factor + 0.5f);
    const int axis_title_font_size = static_cast<int>(base_axis_title_font_size * scale_factor + 0.5f);
    const int label_font_size = static_cast<int>(base_label_font_size * scale_factor + 0.5f);
    const float line_width = base_line_width * scale_factor;
    const float margin = base_margin * scale_factor;

    // Create a single context view - this is the correct approach for VTK charting
    vtkNew<vtkContextView> view;

    // Get the render window and configure it for off-screen rendering
    vtkRenderWindow* renderWindow = view->GetRenderWindow();
    renderWindow->SetOffScreenRendering(1);
    renderWindow->SetSize(width_, height_);

    // Set white background
    view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

    // Create 4 charts
    vtkNew<vtkChartXY> chart_density;
    vtkNew<vtkChartXY> chart_velocity;
    vtkNew<vtkChartXY> chart_pressure;
    vtkNew<vtkChartXY> chart_energy;

    // Add charts to the scene
    view->GetScene()->AddItem(chart_density);
    view->GetScene()->AddItem(chart_velocity);
    view->GetScene()->AddItem(chart_pressure);
    view->GetScene()->AddItem(chart_energy);

    // Calculate chart dimensions
    // Each chart occupies half the width and half the height
    const float w = static_cast<float>(width_);
    const float h = static_cast<float>(height_);
    const float chart_w = w / 2.0f;
    const float chart_h = h / 2.0f;

    // Position charts in 2x2 grid using SetPoint (bottom-left corner) and SetSize
    // VTK scene coordinates: (0,0) is bottom-left

    // Top-left: Density (x: 0 to w/2, y: h/2 to h)
    chart_density->SetAutoSize(false);
    chart_density->SetSize(vtkRectf(margin, chart_h + margin, chart_w - 2 * margin, chart_h - 2 * margin));

    // Top-right: Velocity (x: w/2 to w, y: h/2 to h)
    chart_velocity->SetAutoSize(false);
    chart_velocity->SetSize(vtkRectf(chart_w + margin, chart_h + margin, chart_w - 2 * margin, chart_h - 2 * margin));

    // Bottom-left: Pressure (x: 0 to w/2, y: 0 to h/2)
    chart_pressure->SetAutoSize(false);
    chart_pressure->SetSize(vtkRectf(margin, margin, chart_w - 2 * margin, chart_h - 2 * margin));

    // Bottom-right: Energy (x: w/2 to w, y: 0 to h/2)
    chart_energy->SetAutoSize(false);
    chart_energy->SetSize(vtkRectf(chart_w + margin, margin, chart_w - 2 * margin, chart_h - 2 * margin));

    // Helper lambda to create a data table from DataLayer
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

    // Accessors for each quantity
    auto rho_accessor = [](const DataLayer& dl, int i) { return dl.rho(i); };
    auto u_accessor = [](const DataLayer& dl, int i) { return dl.u(i); };
    auto P_accessor = [](const DataLayer& dl, int i) { return dl.P(i); };
    auto U_accessor = [](const DataLayer& dl, int i) {
        return dl.U(i);
    };  // Specific internal energy

    // Create tables for numerical solution
    auto table_rho = createTable(layer, "Density", rho_accessor);
    auto table_u = createTable(layer, "Velocity", u_accessor);
    auto table_P = createTable(layer, "Pressure", P_accessor);
    auto table_U = createTable(layer, "Energy", U_accessor);

    // Helper lambda to configure chart appearance with scaled fonts
    auto configureChart = [&](vtkChartXY* chart, const std::string& title,
                             const std::string& y_label) {
        chart->SetTitle(title);
        chart->GetTitleProperties()->SetFontSize(title_font_size);
        chart->GetTitleProperties()->SetBold(true);
        chart->GetTitleProperties()->SetColor(0.0, 0.0, 0.0);

        chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("x");
        chart->GetAxis(vtkAxis::BOTTOM)->GetTitleProperties()->SetFontSize(axis_title_font_size);
        chart->GetAxis(vtkAxis::BOTTOM)->GetTitleProperties()->SetColor(0.0, 0.0, 0.0);
        chart->GetAxis(vtkAxis::BOTTOM)->GetLabelProperties()->SetFontSize(label_font_size);
        chart->GetAxis(vtkAxis::BOTTOM)->GetLabelProperties()->SetColor(0.0, 0.0, 0.0);
        chart->GetAxis(vtkAxis::BOTTOM)->GetPen()->SetColor(0, 0, 0);
        chart->GetAxis(vtkAxis::BOTTOM)->GetGridPen()->SetColor(200, 200, 200);

        chart->GetAxis(vtkAxis::LEFT)->SetTitle(y_label);
        chart->GetAxis(vtkAxis::LEFT)->GetTitleProperties()->SetFontSize(axis_title_font_size);
        chart->GetAxis(vtkAxis::LEFT)->GetTitleProperties()->SetColor(0.0, 0.0, 0.0);
        chart->GetAxis(vtkAxis::LEFT)->GetLabelProperties()->SetFontSize(label_font_size);
        chart->GetAxis(vtkAxis::LEFT)->GetLabelProperties()->SetColor(0.0, 0.0, 0.0);
        chart->GetAxis(vtkAxis::LEFT)->GetPen()->SetColor(0, 0, 0);
        chart->GetAxis(vtkAxis::LEFT)->GetGridPen()->SetColor(200, 200, 200);

        chart->SetShowLegend(true);
        // Use integer constants for alignment (2 = RIGHT, 1 = TOP in VTK)
        chart->GetLegend()->SetHorizontalAlignment(2);  // RIGHT
        chart->GetLegend()->SetVerticalAlignment(1);    // TOP
        chart->GetLegend()->GetLabelProperties()->SetFontSize(label_font_size);
    };

    // Helper lambda to add plot to chart with scaled line width
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

    // Configure charts with titles including time
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

    // Add analytical solution first (black line, if available)
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

    // Ensure the render window is properly initialized for off-screen rendering
    renderWindow->SetMultiSamples(0);  // Disable multisampling for off-screen
    renderWindow->SetLineSmoothing(1);
    renderWindow->SetPolygonSmoothing(1);

    // Force the scene to recalculate geometry
    view->GetScene()->SetDirty(true);

    // Render the window (may need multiple passes for VTK charting)
    renderWindow->Render();

    // Second render pass to ensure all charts are drawn
    renderWindow->Render();

    // Capture the image
    vtkNew<vtkWindowToImageFilter> windowToImageFilter;
    windowToImageFilter->SetInput(renderWindow);
    windowToImageFilter->SetScale(1);                // Can increase for higher resolution
    windowToImageFilter->SetInputBufferTypeToRGB();  // Use RGB instead of RGBA
    windowToImageFilter->ReadFrontBufferOff();  // Read from back buffer for off-screen
    windowToImageFilter->Update();

    // Write to PNG file
    std::string filename = GenerateFilename(step);

    vtkNew<vtkPNGWriter> pngWriter;
    pngWriter->SetFileName(filename.c_str());
    pngWriter->SetInputConnection(windowToImageFilter->GetOutputPort());
    pngWriter->Write();
}
