#include "config/ConfigParser.hpp"

#include <algorithm>
#include <stdexcept>

std::optional<bool> ConfigParser::Parse(const std::string& default_config,
                                        const int argc,
                                        char* argv[]) {
    const CliParser cli_parser;
    const std::optional<CliOptions> cli_options =
        cli_parser.Parse(argc, argv, default_config);

    if (!cli_options.has_value()) {
        return std::nullopt;
    }

    config_path_ = cli_options->config_path;
    cli_overrides_ = cli_options->overrides;

    const YamlConfigParser yaml_parser;
    const ParsedYamlConfig parsed = yaml_parser.ParseFile(config_path_);

    defaults_ = parsed.settings;
    initial_conditions_ = parsed.initial_conditions;

    if (cli_options->run_cases.has_value()) {
        run_cases_ = *cli_options->run_cases;
    } else {
        run_cases_ = parsed.run_cases;
    }

    if (run_cases_.empty()) {
        run_cases_.push_back("all");
    }

    return true;
}

const Settings& ConfigParser::GetDefaults() const {
    return defaults_;
}

const std::vector<std::string>& ConfigParser::GetRunCases() const {
    return run_cases_;
}

const std::map<std::string, InitialConditions>& ConfigParser::GetInitialConditions() const {
    return initial_conditions_;
}

const InitialConditions& ConfigParser::GetInitialCondition(const std::string& case_name) const {
    const auto it = initial_conditions_.find(case_name);
    if (it == initial_conditions_.end()) {
        throw std::runtime_error(
                                 "ConfigParser::GetInitialCondition: unknown case '" + case_name + "'"
                                );
    }

    return it->second;
}

Settings ConfigParser::GetCaseSettings(const std::string& case_name) const {
    const auto it = initial_conditions_.find(case_name);
    if (it == initial_conditions_.end()) {
        throw std::runtime_error(
                                 "ConfigParser::GetCaseSettings: unknown case '" + case_name + "'"
                                );
    }

    Settings settings = MergeSettings(defaults_, it->second.overrides, cli_overrides_);
    settings.simulation_case = case_name;
    return settings;
}

std::vector<std::string> ConfigParser::GetAllCaseNames() const {
    std::vector<std::string> names;
    names.reserve(initial_conditions_.size());

    for (const auto& [case_name, _] : initial_conditions_) {
        names.push_back(case_name);
    }

    return names;
}

bool ConfigParser::HasInitialCondition(const std::string& case_name) const {
    return initial_conditions_.find(case_name) != initial_conditions_.end();
}

const std::string& ConfigParser::GetConfigPath() const {
    return config_path_;
}
