{
  description = "CFD Numerical Solver - C++ development environment";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        
        cfd-numerical-solver = pkgs.stdenv.mkDerivation {
          pname = "cfd-numerical-solver";
          version = "0.1.0";
          
          src = ./.;
          
          nativeBuildInputs = with pkgs; [
            cmake
            pkg-config
            doxygen
          ];
          
          buildInputs = with pkgs; [
            xtensor
            xtl
            xsimd
            vtk
            yaml-cpp
          ];
          
          cmakeFlags = [
            "-DCMAKE_BUILD_TYPE=Release"
            "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON"
          ];
          
          # doCheck = true;
          
          meta = with pkgs.lib; {
            description = "CFD Numerical Solver";
            license = licenses.mit;
            platforms = platforms.linux ++ platforms.darwin;
          };
        };
        
      in
      {
        # The buildable package
        packages = {
          default = cfd-numerical-solver;
          cfd-numerical-solver = cfd-numerical-solver;
        };
        
        # Development shell
        devShells.default = pkgs.mkShell {
          inputsFrom = [ cfd-numerical-solver ];
          
          buildInputs = with pkgs; [
            gcc
            cmake
            ninja
            pkg-config
            doxygen
            doxygen_gui
            
            xtensor
            xtl
            xsimd
            vtk
            yaml-cpp
            
            clang-tools
            cppcheck
            gdb
            valgrind
            
            git
          ];
          
          shellHook = ''
            echo "=================================="
            echo "CFD Numerical Solver Dev Environment"
            echo "=================================="
            echo "Compiler: $(gcc --version | head -1)"
            echo "CMake: $(cmake --version | head -1)"
            echo ""
            echo "Quick build commands:"
            echo "  mkdir -p build && cd build"
            echo "  cmake .. && make -j$(nproc)"
            echo "  ./cfd-numerical-solver"
          '';
        };
        
        # Apps for easy execution
        # apps.default = {
        #   type = "app";
        #   program = "${cfd-numerical-solver}/bin/cfd-numerical-solver";
        # };
      });
}
