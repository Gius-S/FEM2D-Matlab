# FEM2D-Laplacian-MATLAB

A MATLAB implementation of the **Finite Element Method (FEM) P1** for solving stationary diffusion problems on 2D non-convex domains.

## Problem

The solver handles the elliptic PDE:

```
-div( κ · ∇u ) = f    in Ω
             u = g    on ∂Ω
```

The included test case solves the **Laplacian on an L-shaped domain** `(-1,1)² \ (-1,0]²`, with exact solution:

```
u(r, θ) = r^(2/3) · sin((2θ + π) / 3)
```

This benchmark is a classic example where the solution has a **corner singularity**, leading to reduced convergence rates:

| Norm | Expected rate |
|------|--------------|
| L²   | O(h^{4/3})  |
| H¹   | O(h^{2/3})  |

## Requirements

- **MATLAB** R2019b or later
- **Partial Differential Equation Toolbox** (for `initmesh` and `decsg`)
- **Symbolic Math Toolbox** (for symbolic differentiation in `Problema.m`)

## Repository Structure

| File | Description |
|------|-------------|
| `Problema.m` | Main script — defines the problem, generates the mesh, runs the solver, and plots convergence results |
| `FEM_2D.m` | Core FEM solver — assembles the stiffness matrix and load vector, applies boundary conditions, solves the linear system, and computes L², H¹, L∞ errors |
| `meshtrans.m` | Converts MATLAB PDE Toolbox mesh format (`P`, `E`, `T`) into the internal struct used by the solver |
| `phih_Pk.m` | Evaluates the reference basis functions φᵢ on the reference triangle for P1, P2, P3 elements |
| `grad_phih_Pk.m` | Evaluates the gradients of the reference basis functions for P1, P2, P3 elements |
| `quadratura.m` | Returns quadrature nodes and weights on the reference triangle for numerical integration |

## How to Run

1. Clone the repository:
   ```bash
   git clone https://github.com/<your-username>/FEM2D-Laplacian-MATLAB.git
   cd FEM2D-Matlab
   ```

2. Open MATLAB and navigate to the project folder.

3. Run the main script:
   ```matlab
   Problema
   ```

This will:
- Build the L-shaped domain and generate a sequence of refined meshes
- Solve the FEM problem on each mesh
- Print the estimated convergence orders in the terminal
- Display two figures:
  - **Figure 1** — log-log convergence plot (L² and H¹ errors vs. mesh size)
  - **Figure 2** — surface plot of the numerical solution on the finest mesh

## Sample Output

```
======== STIMA ORDINI ========
Ordine L2:    1.327     expected: 1.333
Ordine H1:    0.661     expected: 0.667
```

## License

This project is for personal and educational use.
