#ifndef YAMLCONFIGPARSER_HPP
#define YAMLCONFIGPARSER_HPP

#include <map>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "config/InitialConditions.hpp"
#include "config/Settings.hpp"

struct ParsedYamlConfig {
    Settings settings;
    std::vector<std::string> run_cases;
    std::map<std::string, InitialConditions> initial_conditions;
};

class YamlConfigParser {
public:
    YamlConfigParser() = default;

    auto ParseFile(const std::string& filename) -> ParsedYamlConfig;

private:
    static void ParseDefaults(const YAML::Node& defaults_node, Settings& settings);
    static void ParseCases(
        const YAML::Node& cases_node,
        std::map<std::string, InitialConditions>& initial_conditions,
        const Settings& defaults);

    static void ParseRunCases(const YAML::Node& run_node, std::vector<std::string>& run_cases);

    static void ParseMesh(const YAML::Node& node, Settings& settings);
    static void ParsePhysics(const YAML::Node& node, Settings& settings);
    static void ParseChemistry(const YAML::Node& node, Settings& settings);
    static void ParseNumerics(const YAML::Node& node, Settings& settings);
    static void ParseBoundaryConditions(const YAML::Node& node, Settings& settings);
    static void ParseParallel(const YAML::Node& node, Settings& settings);
    static void ParseStopping(const YAML::Node& node, Settings& settings);
    static void ParseLogging(const YAML::Node& node, Settings& settings);
    static void ParseOutput(const YAML::Node& node, Settings& settings);
    static void ParseImmersedBoundaries(const YAML::Node& node, Settings& settings);

    static void ApplyCaseOverrides(const YAML::Node& case_node, InitialConditions& ic);

    static void ParseStructuredInitialCondition(
        const YAML::Node& ic_node,
        InitialConditions& ic,
        const Settings& effective_settings);

    static auto ParseImmersedObjects(const YAML::Node& node) -> std::vector<ImmersedObjectSettings>;

    static auto ReadVectorDouble(const YAML::Node& node) -> std::vector<double>;
    static auto ReadMatrixDouble(const YAML::Node& node) -> std::vector<std::vector<double>>;
    static auto ReadTensorDouble(const YAML::Node& node)
        -> std::vector<std::vector<std::vector<double>>>;

    static void ParseStructured1D(const YAML::Node& ic_node, InitialConditions& ic);
    static void ParseStructured2D(const YAML::Node& ic_node, InitialConditions& ic);
    static void ParseStructured3D(const YAML::Node& ic_node, InitialConditions& ic);

    static void ValidateStructuredShape(const YAML::Node& ic_node, int dim);

    static auto HasKey(const YAML::Node& node, const char* key) -> bool;
};

#endif  // YAMLCONFIGPARSER_HPP