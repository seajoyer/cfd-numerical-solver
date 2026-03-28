#ifndef CLIPARSER_HPP
#define CLIPARSER_HPP

#include <optional>
#include <string>
#include <vector>

#include "config/Settings.hpp"

struct CliOptions {
    std::string config_path;
    std::optional<std::vector<std::string>> run_cases;

    /**
     * Reserved for future CLI overrides.
     */
    CaseSettings overrides;
};

class CliParser {
public:
    CliParser() = default;

    [[nodiscard]] std::optional<CliOptions> Parse(int argc,
                                                  char* argv[],
                                                  const std::string& default_config_path) const;

private:
    [[nodiscard]] static std::vector<std::string> SplitCommaSeparated(const std::string& value);
};

#endif  // CLIPARSER_HPP
