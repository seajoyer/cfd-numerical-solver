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

    [[nodiscard]] ParsedYamlConfig ParseFile(const std::string& filename) const ;

private:
    static constexpr int k_xmin_tag = 0;
    static constexpr int k_xmax_tag = 1;
    static constexpr int k_ymin_tag = 2;
    static constexpr int k_ymax_tag = 3;
    static constexpr int k_zmin_tag = 4;
    static constexpr int k_zmax_tag = 5;

    static void ParseDefaults(const YAML::Node& defaults_node, Settings& settings);
    static void ParseCases(const YAML::Node& cases_node,
                           std::map<std::string, InitialConditions>& initial_conditions,
                           const Settings& defaults);
    static void ParseRunCases(const YAML::Node& run_node, std::vector<std::string>& run_cases);

    static void ParseMesh(const YAML::Node& node, Settings& settings);
    static void ParsePhysics(const YAML::Node& node, Settings& settings);
    static void ParseNumerics(const YAML::Node& node, Settings& settings);
    static void ParseBoundaryConditions(const YAML::Node& node, Settings& settings);
    static void ParseStopping(const YAML::Node& node, Settings& settings);
    static void ParseLogging(const YAML::Node& node, Settings& settings);
    static void ParseOutput(const YAML::Node& node, Settings& settings);
    static void ParseParallel(const YAML::Node& node, Settings& settings);
    static void ParseImmersedBoundaries(const YAML::Node& node, Settings& settings);

    static void ApplyCaseOverrides(const YAML::Node& case_node, InitialConditions& ic);

    static void ParseInitialCondition(const YAML::Node& ic_node,
                                      InitialConditions& ic,
                                      const Settings& effective_settings);

    static void ParseStructuredInitialCondition(const YAML::Node& ic_node,
                                                InitialConditions& ic,
                                                const Settings& effective_settings);

    static void ParseConstantInitialCondition(const YAML::Node& ic_node,
                                              InitialConditions& ic,
                                              const Settings& effective_settings);

    static void ParseStructured1D(const YAML::Node& ic_node, StructuredRegionInitialCondition& ic);
    static void ParseStructured2D(const YAML::Node& ic_node, StructuredRegionInitialCondition& ic);
    static void ParseStructured3D(const YAML::Node& ic_node, StructuredRegionInitialCondition& ic);

    static void ValidateStructuredShape(const YAML::Node& ic_node, int dim);
    static void ValidateSettingsConsistency(const Settings& settings);

    [[nodiscard]] static std::vector<ImmersedObjectSettings> ParseImmersedObjects(const YAML::Node& node);
    [[nodiscard]] static std::vector<double> ReadVectorDouble(const YAML::Node& node);
    [[nodiscard]] static std::vector<std::vector<double>> ReadMatrixDouble(const YAML::Node& node);
    [[nodiscard]] static std::vector<std::vector<std::vector<double>>> ReadTensorDouble(const YAML::Node& node);

    [[nodiscard]] static bool HasKey(const YAML::Node& node, const char* key);
};

#endif  // YAMLCONFIGPARSER_HPP
