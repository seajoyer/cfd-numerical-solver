{
  description = "C++ development environment with";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            libcxx
            libgcc
            cmake

            xtensor
            xtl
            xsimd
            vtk
            yaml-cpp
            doxygen
          ];

          shellHook = ''
            echo "C++ development environment"
            echo "Compiler: $(gcc --version   | head -1)"
            echo "CMake:    $(cmake --version | head -1)"

            # Generate .clangd config for IDE support (clangd from clang-tools)
            cat > .clangd << EOF
CompileFlags:
  Add:
    - "-stdlib=libc++"
    - "-I${pkgs.llvmPackages_latest.libcxx.dev}/include/c++/v1"
  CompilationDatabase: build/
EOF
          '';
        };
      }
    );
}
