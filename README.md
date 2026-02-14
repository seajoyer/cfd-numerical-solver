<div align="center">
    <img src="./doc/full_logo.svg" alt="Project full_logo" width="550" 
        onerror="this.onerror=null; this.src='./full_logo.svg';">
</div>

<br>

<div align="center">
    <pre>Dmitry Sidiuk • Vahe Vahanyan • Kirill Usov</pre>
</div>

<br>

A modular CFD solver for 1D and 2D compressible flow using finite volume methods. Features flexible configuration via YAML, VTK output for visualization in ParaView, and PNG/GIF output for 1D cases. The project is designed for educational and research purposes.

## Features

- **Solvers**: Godunov, Godunov-Kolgan, Godunov-Kolgan-Rodionov (1D), MacCormack (1D), Analytical (1D)
- **Riemann Solvers**: Exact (ideal gas), Acoustic, HLL, HLLC, Osher, Roe, Rusanov
- **Reconstruction**: P0 (piecewise constant), P1 (piecewise linear), ENO, WENO
- **Boundary Conditions**: Free stream, inlet, outlet, reflective, non-reflective, periodic, symmetry, wall
- **Equation of State**: Ideal gas
- **Dimensions**: 1D and 2D
- **Initial Conditions**: Configurable Riemann problems (5 predefined 1D Sod tests), 2D Riemann problems with x-aligned, y-aligned, and four-quadrant discontinuities (Lax-Liu type)
- **Output**: VTK (1D and 2D), PNG and GIF (1D only, with adjustable resolution)

## Dependencies

- C++17 compiler (GCC, Clang)
- CMake 3.25+
- xtensor (auto-fetched if not found)
- yaml-cpp (auto-fetched if not found)
- cxxopts (auto-fetched if not found)
- VTK (auto-fetched if not found)
- Doxygen (optional, auto-fetched if not found)

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

## Configuration

 <details>
 <summary>The <code>config.yaml</code> file defines simulation parameters and initial conditions:</summary>

```yaml
config:
    run_cases: [sod1, sod2]  # or simply: all

    global:
        solver: godunov
        riemann_solver: exact
        reconstruction: P0
        time_integrator: euler
        left_boundary: free_stream
        right_boundary: free_stream
        bottom_boundary: free_stream   # 2D: y-min boundary
        top_boundary: free_stream      # 2D: y-max boundary

        global_limiter: false
        diffusion: false
        viscosity: false

        N: 1000
        Nx: 0      # 2D: cells in x (0 = use N)
        Ny: 0      # 2D: cells in y (0 = use N)
        cfl: 0.5
        padding: 2
        gamma: 1.4
        dim: 1
        L_x: 1
        L_y: 1
        L_z: 1

        Q_user: 2

        x0: 0.5
        y0: 0.5    # 2D: discontinuity y-position
        analytical: true

        t_end: 0.25
        step_end: 0

        log_every_steps: 50
        log_every_time: 0.0

        output_every_steps: 50
        output_every_time: 0.0

        output_formats: [vtk, png, gif2560x1600]
        output_dir: "../result"

    # Each case can override any global setting
    cases:
        sod1:
            rho_L: 1.0
            u_L: 0.0
            P_L: 1.0
            rho_R: 0.125
            u_R: 0.0
            P_R: 0.1
            x0: 0.5
            t_end: 0.25

        # 2D four-quadrant Riemann problem (Lax-Liu Config 3)
        lax_liu_3:
            dim: 2
            N: 400
            ic_type: quadrant
            x0: 0.5
            y0: 0.5
            cfl: 0.4
            t_end: 0.3
            rho_Q1: 1.5
            u_Q1: 0.0
            v_Q1: 0.0
            P_Q1: 1.5
            rho_Q2: 0.5323
            u_Q2: 1.206
            v_Q2: 0.0
            P_Q2: 0.3
            rho_Q3: 0.138
            u_Q3: 1.206
            v_Q3: 1.206
            P_Q3: 0.029
            rho_Q4: 0.5323
            u_Q4: 0.0
            v_Q4: 1.206
            P_Q4: 0.3
            output_formats: [vtk]
            output_every_steps: 200
```
</details>

Override configuration parameters as follows:
```bash
./cfd-numerical-solver --run-cases sod1,sod3 --N-cells 2000 --cfl 0.4
```

Check `./cfd-numerical-solver --help` for a full list of command-line options, including flags to list supported solvers, Riemann solvers, boundary conditions, simulation cases, and output formats.

## Usage

Run with default configuration:
```bash
./cfd-numerical-solver
```

The solver reads `../config.yaml` by default and outputs to `../result/`.

### 1D Examples

```bash
# Classic Sod shock tube with analytical comparison
./cfd-numerical-solver --run-cases sod1 --N-cells 1000 --analytical true

# Strong shock (Sod 3) with HLLC Riemann solver
./cfd-numerical-solver --run-cases sod3 --riemann-solver hllc --cfl 0.4
```

### 2D Examples

```bash
# 2D Sod shock tube along x-axis
./cfd-numerical-solver --run-cases sod1_2d --dim 2 --N-cells 200

# Lax-Liu quadrant problem
./cfd-numerical-solver --run-cases lax_liu_3 --N-cells 400 --cfl 0.4
```

## 2D Solver Details

The 2D extension uses an unsplit finite volume method with dimensionally-independent flux sweeps:

- **Spatial discretization**: Godunov's method with P0 (piecewise constant) reconstruction. The x-sweep and y-sweep each reuse the existing 1D Riemann solvers via velocity rotation: the y-sweep passes `v` as the normal velocity and passively advects `u`.
- **Time integration**: Forward Euler (explicit). Higher-order time integrators (SSPRK2, SSPRK3) are available for 1D only.
- **CFL condition**: Uses the unsplit 2D spectral radius: `dt = CFL / ((|u|+c)/dx + (|v|+c)/dy)`.
- **Boundary conditions**: All 8 boundary types support 2D with correct axis-dependent normal/tangential velocity handling.
- **Initial conditions**: Three IC types for 2D — `x_riemann` (vertical interface), `y_riemann` (horizontal interface), and `quadrant` (four-quadrant Lax-Liu type).
- **Output**: VTK structured points format, viewable in ParaView. PNG/GIF output is not yet supported for 2D.

### Current 2D Limitations

- Only P0 (piecewise constant) reconstruction is supported; higher-order spatial reconstruction is planned.
- Time integration is limited to Forward Euler; SSPRK2/SSPRK3 are not yet available for 2D.
- MacCormack and Godunov-Kolgan-Rodionov solvers fall back to first-order Godunov in 2D.
- Analytical solution comparison is not available for 2D cases.
- Diffusion filter and global limiter are not applied in 2D.
- Visual output (PNG, GIF) is 1D-only; use VTK + ParaView for 2D visualization.

## Output Structure

Results are organized hierarchically:
```
result/
├── run_01-12-2025_10:29:00:018
│   ├── sod1
│   │   ├── godunov__R_p0__N_1000__CFL_5e-1.gif
│   │   ├── png
│   │   │   ├── step_0000.png
│   │   │   ├── step_0050.png
│   │   │   └── ...
│   │   └── vtk
│   │       ├── analytical
│   │       │   ├── step_0000.vtk
│   │       │   └── ...
│   │       └── godunov__R_p0__N_1000__CFL_5e-1
│   │           ├── godunov__R_p0__N_1000__CFL_5e-1__step_0000.vtk
│   │           └── ...
│   ├── lax_liu_3
│   │   └── vtk
│   │       └── godunov__R_p0__N_400x400__CFL_4e-1
│   │           ├── godunov__R_p0__N_400x400__CFL_4e-1__step_0000.vtk
│   │           └── ...
│   └── ...
│
└── run_16-12-2025_16:57:02:437
    └── ...
```

## License

Licensed under the Apache License, Version 2.0. See http://www.apache.org/licenses/LICENSE-2.0 for details.
