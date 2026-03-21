#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkRenderingContextOpenGL2);

#include <exception>
#include <iostream>

#include "RunManager.hpp"
#include "parallel/MPIContext.hpp"

namespace {
    class MPIInitializerGuard {
    public:
        MPIInitializerGuard(int& argc, char**& argv) {
            MPIContext::Initialize(argc, argv);
        }

        ~MPIInitializerGuard() {
            MPIContext::Finalize();
        }

        MPIInitializerGuard(const MPIInitializerGuard&) = delete;
        auto operator=(const MPIInitializerGuard&) -> MPIInitializerGuard& = delete;
    };
}  // namespace

auto main(int argc, char* argv[]) -> int {
    MPIInitializerGuard mpi_guard(argc, argv);

    try {
        RunManager run_manager;
        return run_manager.Run(argc, argv);
    }
    catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << '\n';
        return 1;
    }
    catch (...) {
        std::cerr << "Fatal error: unknown exception\n";
        return 1;
    }
}