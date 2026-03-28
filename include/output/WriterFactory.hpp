#ifndef WRITERFACTORY_HPP
#define WRITERFACTORY_HPP

#include <memory>
#include <string>
#include <vector>

class StepWriter;

/**
 * @class WriterFactory
 * @brief Factory for creating output writers.
 *
 * Currently supported formats:
 *  - "vtk"
 */
class WriterFactory final {
public:
    /**
     * @brief Create one writer for the specified output format.
     * @param output_format Output format identifier.
     * @param output_dir Directory where files will be written.
     * @return Unique pointer to created writer.
     */
    [[nodiscard]] static std::unique_ptr<StepWriter> Create(const std::string& output_format,
                                                            const std::string& output_dir);

    /**
     * @brief Create multiple writers for the specified output formats.
     * @param output_formats List of output format identifiers.
     * @param output_dir Directory where files will be written.
     * @return Vector of created writers.
     */
    [[nodiscard]] static std::vector<std::unique_ptr<StepWriter>> CreateMultiple(
        const std::vector<std::string>& output_formats,
        const std::string& output_dir
    );

    /**
     * @brief Check whether output format is supported.
     */
    [[nodiscard]] static bool IsFormatSupported(const std::string& format);

    /**
     * @brief Get list of supported output formats.
     */
    [[nodiscard]] static std::vector<std::string> GetSupportedFormats();
};

#endif  // WRITERFACTORY_HPP
