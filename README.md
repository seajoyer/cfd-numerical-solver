<div align="center">
    <img src="./doc/full_logo.svg" alt="Project full_logo" width="550" 
        onerror="this.onerror=null; this.src='./full_logo.svg';">
</div>

<br>

<div align="center">
    <pre>Dmitry Sidiuk • Vahe Vahanyan • Kirill Usov</pre>
</div>

<br>

A modular CFD solver for 1D compressible flow using finite volume methods. Features flexible configuration via YAML and VTK output for visualization in ParaView. The project is designed for educational purposes.

## Features

- **Solvers**: Godunov, Godunov-Kolgan, Godunov-Kolgan-Rodionov, Analytical
- **Riemann Solvers**: Exact (ideal gas), Acoustic, HLL, HLLC
- **Reconstruction**: P0 (piecewise constant), P1 (piecewise linear), ENO, WENO
- **Boundary Conditions**: Free stream, inlet, outlet, reflective, non-reflective, periodic, symmetry, wall
- **Equation of State**: Ideal gas
- **Initial Conditions**: Configurable Riemann problems (5 predefined Sod tests + custom cases)
- **Dimensions**: Currently 1D (2D/3D planned)
- **Output**: VTK files with organized directory structure

## Dependencies

- C++17 compiler (GCC, Clang)
- CMake 3.25+
- xtensor (auto-fetched if not found)
- yaml-cpp (auto-fetched if not found)
- cxxopts (auto-fetched if not found)
- VTK (required, install via package manager)
- Doxygen (optional, for documentation)

## Building from Source

1. Clone the repository:
   ```bash
   git clone https://github.com/seajoyer/cfd-numerical-solver
   cd cfd-numerical-solver
   ```

2. Create build directory:
   ```bash
   mkdir build && cd build
   ```

3. Configure and build:
   ```bash
   cmake ..
   cmake --build .
   ```

4. (Optional) Generate documentation:
   ```bash
   cmake --build . --target docs
   ```
   Documentation available at `build/docs/html/index.html`

## Usage

Run with default configuration:
```bash
./cfd-numerical-solver
```

The solver reads `../config.yaml` by default and outputs to `../data/`.

### Command Line Options

Override configuration parameters:
```bash
./cfd-numerical-solver -i sod3 --N-cells 2000 --cfl 0.4
```

Key options:
- `-i, --simulation-case`: Select simulation case (e.g., `sod3`, `custom_case`, or `all`)
- `-N, --N-cells`: Number of grid cells
- `--cfl`: CFL number
- `-s, --solver`: Solver type
- `--riemann-solver`: Riemann solver type
- `--reconstruction`: Reconstruction scheme
- `-o, --output-dir`: Output directory
- `-h, --help`: Show all options

## Configuration

 <details>
 <summary>The <code>config.yaml</code> file defines simulation parameters and initial conditions:</summary>
```yaml
config:
  settings:
    solver: godunov
    riemann_solver: exact
    reconstruction: P0
    left_boundary: free_stream
    right_boundary: free_stream

    N: 1000
    cfl: 0.5
    padding: 2
    gamma: 1.4
    dim: 1
    L_x: 1

    Q_user: 2

    t_end: 0.25
    step_end: 0

    log_every_steps: 50
    log_every_time: 0.0

    output_every_steps: 50
    output_every_time: 0.0

    output_format: vtk
    output_dir: "../data"

    # Use "all" to run all available cases
    simulation_case: sod1
    x0: 0.5
    analytical: true

  initial_conditions:
    sod1:
      rho_L: 1.0
      u_L: 0.0
      P_L: 1.0
      rho_R: 0.125
      u_R: 0.0
      P_R: 0.1
      x0: 0.5
      t_end: 0.25
    
    sod2:
      rho_L: 1.0
      u_L: -2.0
      P_L: 0.4
      rho_R: 1.0
      u_R: 2.0
      P_R: 0.4
      x0: 0.5
      t_end: 0.15
    
    # Add custom cases here
    # custom_case:
    #   rho_L: 2.0
    #   u_L: 1.0
    #   P_L: 5.0
    #   rho_R: 0.5
    #   u_R: -0.5
    #   P_R: 2.0
    #   x0: 0.3
    #   t_end: 0.1
```

</details>

### Output Structure

Results are organized hierarchically:

```
data/
├── sod1/  # сase names are added only if `simulation-case` is `all`.
│   ├── godunov__R_p0__N_1000__CFL_5e-1/
│   │   └── *.vtk
│   └── analytical/
│       └── *.vtk
├── sod3/
│   └── ...
```

Each numerical simulation creates a directory with solver parameters in the name, while analytical solutions go to `analytical/` subdirectories.

### Running Multiple Cases

Set `simulation_case: all` to run all defined initial conditions in sequence:

```bash
./cfd-numerical-solver -i all
```

## License

Licensed under the Apache License, Version 2.0. See http://www.apache.org/licenses/LICENSE-2.0 for details.
