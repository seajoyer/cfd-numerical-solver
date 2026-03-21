#ifndef CLIPARSER_HPP
#define CLIPARSER_HPP

#include <optional>
#include <string>
#include <vector>

struct CliOptions {
    std::string config_path;
    std::optional<std::vector<std::string>> run_cases;
};

class CliParser {
public:
    CliParser() = default;

    auto Parse(int argc, char* argv[], const std::string& default_config_path)
        -> std::optional<CliOptions>;

private:
    static auto SplitCommaSeparated(const std::string& value) -> std::vector<std::string>;
};

#endif  // CLIPARSER_HPP