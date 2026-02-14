<div align="center">
    <img src="./doc/full_logo.svg" alt="Project full_logo" width="550" 
        onerror="this.onerror=null; this.src='./full_logo.svg';">
</div>

<br>

<div align="center">
    <pre>Dmitry Sidiuk • Vahe Vahanyan • Kirill Usov</pre>
</div>

<br>

A modular CFD solver for 1D compressible flow using finite volume methods. Features flexible configuration via YAML with support of VTK output for visualization in ParaView. The project is designed for educational purposes.

## Features

- **Solvers**: Godunov, Godunov-Kolgan, Godunov-Kolgan-Rodionov, McCormak, Analytical
- **Riemann Solvers**: Exact (ideal gas), Acoustic, HLL, HLLC, Osher, Roe, Rusanov
- **Reconstruction**: P0 (piecewise constant), P1 (piecewise linear), ENO, WENO
- **Boundary Conditions**: Free stream, inlet, outlet, reflective, non-reflective, periodic, symmetry, wall
- **Equation of State**: Ideal gas
- **Initial Conditions**: Configurable Riemann problems (5 predefined Sod tests + custom cases)
- **Dimensions**: Currently 1D (2D/3D planned)
- **Output**: VTK, PNG and GIF (with adjustable resolution)

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

        global_limiter: false
        diffusion: true
        viscosity: false

        N: 1000
        cfl: 0.5
        padding: 2
        gamma: 1.4
        dim: 1
        L_x: 1
        L_y: 1
        L_z: 1

        Q_user: 2

        x0: 0.5
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
        
        sod2:
            rho_L: 1.0
            u_L: -2.0
            P_L: 0.4
            rho_R: 1.0
            u_R: 2.0
            P_R: 0.4

            x0: 0.5
            t_end: 0.15
        
        sod3:
            rho_L: 1.0
            u_L: 0.0
            P_L: 1000.0
            rho_R: 1.0
            u_R: 0.0
            P_R: 0.01

            x0: 0.5
            t_end: 0.012
        
        sod4:
            rho_L: 1.0
            u_L: 0.0
            P_L: 0.01
            rho_R: 1.0
            u_R: 0.0
            P_R: 100.0

            x0: 0.5
            t_end: 0.035
        
        sod5:
            rho_L: 5.99924
            u_L: 19.5975
            P_L: 460.894
            rho_R: 5.99242
            u_R: -6.19633
            P_R: 46.0950

            x0: 0.5
            t_end: 0.035
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

## Output Structure

Results are organized hierarchically:

```
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
│   │       │   ├── step_0050.vtk
│   │       │   └── ...
│   │       └── godunov__R_p0__N_1000__CFL_5e-1
│   │           ├── godunov__R_p0__N_1000__CFL_5e-1__step_0000.vtk
│   │           ├── godunov__R_p0__N_1000__CFL_5e-1__step_0050.vtk
│   │           └── ...
│   ├── custom-case
│   │   └── ...
│   └── ...
│
└── run_16-12-2025_16:57:02:437
    └── ...
```

## License

Licensed under the Apache License, Version 2.0. See http://www.apache.org/licenses/LICENSE-2.0 for details.
