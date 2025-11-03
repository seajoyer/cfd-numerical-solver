#ifndef VTKWRITER_HPP
#define VTKWRITER_HPP

#include "../config/StepWriter.hpp"
#include <string>

class VtkWriter : public StepWriter {
public:
    explicit VtkWriter(std::string outputDir);

    void Write(const DataLayer &layer, std::size_t step, double time) const override;

private:
    std::string outputDir;
};

#endif // VTKWRITER_HPP