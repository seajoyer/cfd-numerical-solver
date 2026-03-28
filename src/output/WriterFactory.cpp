#include "output/WriterFactory.hpp"

#include <algorithm>
#include <cctype>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "output/StepWriter.hpp"
#include "output/VTKWriter.hpp"

namespace {
    [[nodiscard]] std::string ToLower(std::string value) {
        std::transform(value.begin(), value.end(), value.begin(),
                       [](const unsigned char c) {
                           return static_cast<char>(std::tolower(c));
                       });
        return value;
    }
} // namespace

std::unique_ptr<StepWriter> WriterFactory::Create(const std::string& output_format,
                                                  const std::string& output_dir) {
    const std::string format_lower = ToLower(output_format);

    if (format_lower == "vtk") {
        return std::make_unique<VTKWriter>(output_dir);
    }

    throw std::runtime_error(
                             "WriterFactory::Create: unknown output format: " + output_format +
                             ". Supported formats: vtk"
                            );
}

std::vector<std::unique_ptr<StepWriter>> WriterFactory::CreateMultiple(
    const std::vector<std::string>& output_formats,
    const std::string& output_dir
) {
    std::vector<std::unique_ptr<StepWriter>> writers;
    writers.reserve(output_formats.size());

    for (const std::string& format : output_formats) {
        writers.push_back(Create(format, output_dir));
    }

    return writers;
}

bool WriterFactory::IsFormatSupported(const std::string& format) {
    return ToLower(format) == "vtk";
}

std::vector<std::string> WriterFactory::GetSupportedFormats() {
    return {"vtk"};
}
