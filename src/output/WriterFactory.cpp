#include "output/WriterFactory.hpp"

#include <algorithm>
#include <stdexcept>

#include "output/GIFWriter.hpp"
#include "output/PNGWriter.hpp"
#include "output/VTKWriter.hpp"

void WriterFactory::ParseResolution(const std::string& format_lower,
                                    const std::string& prefix,
                                    int& width,
                                    int& height) {
    // Default resolution
    width = 1200;
    height = 900;

    // Check if resolution is specified after the prefix
    if (format_lower.length() > prefix.length()) {
        std::size_t x_pos = format_lower.find('x', prefix.length());
        if (x_pos != std::string::npos) {
            try {
                width = std::stoi(format_lower.substr(prefix.length(), x_pos - prefix.length()));
                height = std::stoi(format_lower.substr(x_pos + 1));

                // Validate resolution bounds
                if (width < 100 || width > 10000 || height < 100 || height > 10000) {
                    throw std::runtime_error(
                        prefix + " resolution must be between 100x100 and 10000x10000");
                }
            } catch (const std::invalid_argument&) {
                throw std::runtime_error(
                    "Invalid " + prefix + " resolution format. Use: " + prefix +
                    "<width>x<height> (e.g., " + prefix + "1920x1080)");
            } catch (const std::out_of_range&) {
                throw std::runtime_error(prefix + " resolution values out of range");
            }
        }
    }
}

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

    // Parse PNG format: "png" or "png<width>x<height>"
    if (format_lower.substr(0, 3) == "png") {
        int width = 0;
        int height = 0;
        ParseResolution(format_lower, "png", width, height);

        // PNG writer doesn't need is_analytical flag - it handles both in same file
        return std::make_unique<PNGWriter>(output_dir, width, height);
    }

    // Parse GIF format: "gif" or "gif<width>x<height>"
    if (format_lower.substr(0, 3) == "gif") {
        int width = 0;
        int height = 0;
        ParseResolution(format_lower, "gif", width, height);

        // GIF writer stores to the case directory directly (not a subdirectory)
        // The output_dir passed here should be the case directory
        return std::make_unique<GIFWriter>(output_dir, width, height);
    }

    throw std::runtime_error("Unknown output format: " + output_format +
                             ". Supported formats: vtk, png[<width>x<height>], "
                             "gif[<width>x<height>]");
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

    if (format_lower == "vtk") {
        return true;
    }

    if (format_lower.starts_with("png")) {
        return true;
    }

    if (format_lower.starts_with("gif")) {
        return true;
    }

    return false;
}

auto WriterFactory::GetSupportedFormats() -> std::vector<std::string> {
    return {"vtk", "png", "png<width>x<height>", "gif", "gif<width>x<height>"};
}
