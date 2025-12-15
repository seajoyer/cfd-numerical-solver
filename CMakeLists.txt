cmake_minimum_required(VERSION 3.25)

set(NAME cfd-numerical-solver)
project(${NAME} LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(FETCHCONTENT_QUIET ON)
set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE ALWAYS)

include(FetchContent)

include(cmake/Dependencies.cmake)
include(cmake/Documentation.cmake)

# Explicitly list your sources
add_executable(${NAME}
    src/main.cpp

    src/Simulation.cpp

    src/solver/SolverFactory.cpp
    src/solver/GodunovSolver.cpp
    src/solver/GodunovKolganRodionovSolver.cpp
    src/solver/FiniteVolumeSolver.cpp
    src/solver/MacCormackSolver.cpp
    src/solver/AnalyticalSolver.cpp
    include/solver/Solver.hpp
    src/solver/EOS.cpp
    src/solver/TimeStepCalculator.cpp
    src/solver/PositivityLimiter.cpp

    src/time/ForwardEulerTimeIntegrator.cpp
    src/time/MacCormackTimeIntegrator.cpp
    src/time/SSPRK3TimeIntegrator.cpp
    src/time/SSPRK2TimeIntegrator.cpp
    src/time/TimeIntegrator.cpp

    src/limiter/GlobalLimiter.cpp
    src/limiter/VacuumFixLimiter.cpp

    src/filter/SolutionFilter.cpp

    src/viscosity/VNRArtificialViscosity.cpp
    include/viscosity/ArtificialViscosity.hpp

    src/spatial/ForwardEulerSpatialOperator.cpp
    src/spatial/BackwardEulerSpatialOperator.cpp
    src/spatial/GodunovSpatialOperator.cpp
    src/spatial/GodunovSpatialOperator.cpp
    src/spatial/GodunocKolganRodionovSpatialOperator.cpp
    include/spatial/SpatialOperator.hpp

    src/reconstruction/P0Reconstruction.cpp
    src/reconstruction/P1Reconstruction.cpp
    src/reconstruction/ENOReconstruction.cpp
    src/reconstruction/WENOReconstruction.cpp
    include/reconstruction/Reconstruction.hpp

    src/riemann/RusanovRiemannSolver.cpp
    src/riemann/HLLRiemannSolver.cpp
    src/riemann/HLLCRiemannSolver.cpp
    src/riemann/RoeRiemannSolver.cpp
    src/riemann/OsherRiemannSolver.cpp
    src/riemann/ExactIdealGasRiemannSolver.cpp
    src/riemann/AcousticRiemannSolver.cpp

    src/output/WriterFactory.cpp
    src/output/VTKWriter.cpp
    src/output/PNGWriter.cpp
    src/output/GIFWriter.cpp

    src/data/DataLayer.cpp
    src/data/Variables.cpp

    src/config/ConfigParser.cpp
    include/config/Settings.hpp
    include/config/InitialConditions.hpp

    src/bc/BoundaryFactory.cpp
    src/bc/BoundaryManager.cpp
    src/bc/FreeStreamBoundary.cpp
    src/bc/InletBoundary.cpp
    src/bc/OutletBoundary.cpp
    src/bc/PeriodicBoundary.cpp
    src/bc/ReflectiveBoundary.cpp
    src/bc/NonReflectiveBoundary.cpp
    src/bc/SymmetryBoundary.cpp
    src/bc/WallBoundary.cpp
    include/bc/BoundaryCondition.hpp

    src/utils/StringUtils.cpp
)

target_include_directories(${NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/include/external
)

target_link_libraries(${NAME}
    PRIVATE
        xtensor
        yaml-cpp::yaml-cpp
        cxxopts::cxxopts

        VTK::CommonCore
        VTK::CommonDataModel
        VTK::IOLegacy
        VTK::ChartsCore
        VTK::ViewsContext2D
        VTK::RenderingContext2D
        VTK::IOImage
        VTK::RenderingContextOpenGL2
        VTK::RenderingOpenGL2
)
