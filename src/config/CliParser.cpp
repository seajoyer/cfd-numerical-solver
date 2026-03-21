#include "config/CliParser.hpp"

#include <cxxopts.hpp>
#include <iostream>
#include <sstream>

auto CliParser::Parse(int argc, char* argv[], const std::string& default_config_path)
    -> std::optional<CliOptions> {
    try {
        cxxopts::Options options(
            "cfd_numerical_solver",
            "CFD Numerical Solver"
        );

        options.add_options()
            ("h,help", "Print help")
            ("c,config", "Path to config file", cxxopts::value<std::string>())
            ("run-cases", "Comma-separated list of cases to run", cxxopts::value<std::string>());

        const auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help() << '\n';
            return std::nullopt;
        }

        CliOptions cli;
        cli.config_path = default_config_path;

        if (result.count("config")) {
            cli.config_path = result["config"].as<std::string>();
        }

        if (result.count("run-cases")) {
            cli.run_cases = SplitCommaSeparated(result["run-cases"].as<std::string>());
        }

        return cli;
    }
    catch (const cxxopts::exceptions::exception& e) {
        std::cerr << "CLI parsing error: " << e.what() << '\n';
        return std::optional<CliOptions>{};
    }
}

auto CliParser::SplitCommaSeparated(const std::string& value) -> std::vector<std::string> {
    std::vector<std::string> result;
    std::istringstream ss(value);
    std::string item;

    while (std::getline(ss, item, ',')) {
        const auto begin = item.find_first_not_of(" \t");
        if (begin == std::string::npos) {
            continue;
        }

        const auto end = item.find_last_not_of(" \t");
        result.push_back(item.substr(begin, end - begin + 1));
    }

    return result;
}