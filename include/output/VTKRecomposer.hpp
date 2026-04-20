#ifndef VTKRECOMPOSER_HPP
#define VTKRECOMPOSER_HPP

#include <string>
#include <vector>

class VTKRecomposer final {
public:
    static void RecomposeCaseDirectory(const std::string& case_dir,
                                       int expected_rank_count);

private:
    [[nodiscard]] static std::vector<std::string> FindRankDirectories(const std::string& case_dir,
                                                                      int expected_rank_count);
};

#endif  // VTKRECOMPOSER_HPP
