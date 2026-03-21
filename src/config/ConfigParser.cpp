#include "config/ConfigParser.hpp"

#include <algorithm>
#include <stdexcept>

auto ConfigParser::Parse(const std::string& default_config, int argc, char* argv[])
    -> std::optional<bool> {
    CliParser cli_parser;
    const auto cli_result = cli_parser.Parse(argc, argv, default_config);

    if (!cli_result.has_value()) {
        return std::nullopt;
    }

    const CliOptions cli = *cli_result;
    config_path_ = cli.config_path;

    try {
        YamlConfigParser yaml_parser;
        ParsedYamlConfig parsed = yaml_parser.ParseFile(config_path_);

        settings_ = parsed.settings;
        initial_conditions_ = std::move(parsed.initial_conditions);
        run_cases_ = std::move(parsed.run_cases);

        if (cli.run_cases.has_value()) {
            run_cases_ = *cli.run_cases;
        }

        return true;
    }
    catch (const std::exception&) {
        throw;
    }
}

auto ConfigParser::GetSettings() const -> const Settings& {
    return settings_;
}

auto ConfigParser::GetRunCases() const -> const std::vector<std::string>& {
    return run_cases_;
}

auto ConfigParser::GetInitialConditions() const
    -> const std::map<std::string, InitialConditions>& {
    return initial_conditions_;
}

auto ConfigParser::GetInitialCondition(const std::string& case_name) const
    -> const InitialConditions& {
    const auto it = initial_conditions_.find(case_name);
    if (it == initial_conditions_.end()) {
        throw std::out_of_range("Case '" + case_name + "' not found");
    }
    return it->second;
}

auto ConfigParser::GetCaseSettings(const std::string& case_name) const -> Settings {
    const InitialConditions& ic = GetInitialCondition(case_name);
    return MergeSettings(settings_, ic.overrides, CaseSettings{});
}

auto ConfigParser::GetAllCaseNames() const -> std::vector<std::string> {
    std::vector<std::string> names;
    names.reserve(initial_conditions_.size());

    for (const auto& [name, _] : initial_conditions_) {
        names.push_back(name);
    }

    std::sort(names.begin(), names.end());
    return names;
}

auto ConfigParser::HasInitialCondition(const std::string& case_name) const -> bool {
    return initial_conditions_.find(case_name) != initial_conditions_.end();
}

auto ConfigParser::GetConfigPath() const -> const std::string& {
    return config_path_;
}