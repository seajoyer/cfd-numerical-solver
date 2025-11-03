#include "visualization/VTKWriter.hpp"
#include "data/DataLayer.hpp"
#include <utility>

VtkWriter::VtkWriter(std::string outputDir)
    : outputDir(std::move(outputDir)) {
}

void VtkWriter::Write(const DataLayer & /*layer*/, std::size_t /*step*/, double /*time*/) const {
    // TODO: запись данных шага в VTK/VTU файл (пока заглушка)
}