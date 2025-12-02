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

    if (format_lower == "png") {
        // PNG writer doesn't need is_analytical flag - it handles both in same file
        return std::make_unique<PNGWriter>(output_dir);
    }

    // Future formats can be added here:
    //
    // if (format_lower == "plt") {
    //     return std::make_unique<PLTWriter>(output_dir);
    // }

    throw std::runtime_error("Unknown output format: " + output_format + 
                            ". Supported formats: vtk, png");
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
    
    return (format_lower == "vtk" || format_lower == "png");
}

auto WriterFactory::GetSupportedFormats() -> std::vector<std::string> {
    return {"vtk", "png"};
}
