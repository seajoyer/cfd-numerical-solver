#include "output/WriterFactory.hpp"

#include <algorithm>
#include <stdexcept>

#include "output/PNGWriter.hpp"
#include "output/VTKWriter.hpp"

auto WriterFactory::Create(const std::string& output_format,
                           const std::string& output_dir,
                           bool is_analytical) -> std::unique_ptr<StepWriter> {
    // Convert to lowercase for case-insensitive comparison
    std::string format_lower = output_format;
    std::transform(format_lower.begin(), format_lower.end(), format_lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (format_lower == "vtk") {
        return std::make_unique<VTKWriter>(output_dir, is_analytical);
    }

    // Parse PNG format: "png" or "png<width>x<height>" (e.g., "png1200x900")
    if (format_lower.substr(0, 3) == "png") {
        int width = 1200;   // Default width
        int height = 900;   // Default height
        
        // Parse resolution if specified: png1200x900
        if (format_lower.length() > 3) {
            std::size_t x_pos = format_lower.find('x', 3);
            if (x_pos != std::string::npos) {
                try {
                    width = std::stoi(format_lower.substr(3, x_pos - 3));
                    height = std::stoi(format_lower.substr(x_pos + 1));
                    
                    // Validate resolution bounds
                    if (width < 100 || width > 10000 || height < 100 || height > 10000) {
                        throw std::runtime_error("PNG resolution must be between 100x100 and 10000x10000");
                    }
                } catch (const std::invalid_argument&) {
                    throw std::runtime_error("Invalid PNG resolution format. Use: png<width>x<height> (e.g., png1920x1080)");
                } catch (const std::out_of_range&) {
                    throw std::runtime_error("PNG resolution values out of range");
                }
            }
        }
        
        // PNG writer doesn't need is_analytical flag - it handles both in same file
        return std::make_unique<PNGWriter>(output_dir, width, height);
    }

    // Future formats can be added here:
    //
    // if (format_lower == "plt") {
    //     return std::make_unique<PLTWriter>(output_dir);
    // }

    throw std::runtime_error("Unknown output format: " + output_format + 
                            ". Supported formats: vtk, png[<width>x<height>]");
}

auto WriterFactory::CreateMultiple(const std::vector<std::string>& output_formats,
                                   const std::string& output_dir,
                                   bool is_analytical)
    -> std::vector<std::unique_ptr<StepWriter>> {
    std::vector<std::unique_ptr<StepWriter>> writers;
    writers.reserve(output_formats.size());

    for (const auto& format : output_formats) {
        writers.push_back(Create(format, output_dir, is_analytical));
    }

    return writers;
}

auto WriterFactory::IsFormatSupported(const std::string& format) -> bool {
    std::string format_lower = format;
    std::transform(format_lower.begin(), format_lower.end(), format_lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    
    // Support "vtk" or "png" or "png<width>x<height>"
    if (format_lower == "vtk") {
        return true;
    }
    
    if (format_lower.substr(0, 3) == "png") {
        return true;
    }
    
    return false;
}

auto WriterFactory::GetSupportedFormats() -> std::vector<std::string> {
    return {"vtk", "png", "png<width>x<height>"};
}
