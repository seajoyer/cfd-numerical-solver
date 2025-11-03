{
  description = "C++ development environment with";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            # Compiler and build tools
            gcc
            cmake
            ninja
            pkg-config

            # xtensor and dependencies
            xtensor
            xtl
            xsimd
            vtk

            # Optional: additional useful tools
            # clang-tools
            cppcheck
            gdb
            valgrind
          ];

          shellHook = ''
            echo "C++ development environment with xtensor"
            echo "Compiler: $(gcc --version | head -1)"
            echo "CMake: $(cmake --version | head -1)"
          '';
        };

        # Optional: You can also define a package if you want to build your project
        packages.default = pkgs.stdenv.mkDerivation {
          name = "xtensor-project";
          src = self;
          nativeBuildInputs = with pkgs; [ cmake ninja pkg-config ];
          buildInputs = with pkgs; [ xtensor xtl xsimd ];
          cmakeFlags = [ "-DCMAKE_BUILD_TYPE=Release" ];
        };
      });
}
