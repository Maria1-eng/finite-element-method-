# Finite Element Method

Two implementations of the Finite Element Method (FEM) applied to different physical problems, developed as part of a Finite Elements course at Universidad Tecnológica de la Mixteca.

---

## 1. Poisson Equation on a Square Plate — Python

### Problem

Solve the 2D Poisson equation on a square domain with homogeneous Dirichlet boundary conditions:

$$-\left(\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}\right) = f_0, \qquad u = 0 \text{ on } \Gamma$$

The diagonal symmetry of the domain (x = y) allows reducing the problem to a triangular region discretized with **15 nodes** and **15 linear triangular elements**.

### Rayleigh-Ritz Method

Multiplying by a test function φ and integrating over the domain, then applying the divergence theorem:

$$\iint \left[\frac{\partial\phi}{\partial x}\frac{\partial u}{\partial x} + \frac{\partial\phi}{\partial y}\frac{\partial u}{\partial y}\right]dxdy = \iint \phi\, f_0\, dxdy$$

Using linear shape functions N₁, N₂, N₃ over each triangular element and evaluating the integrals yields the element stiffness matrix.

### Master Element

For the triangular elements in this mesh (with f₀ = 1):

$$\frac{1}{2}\begin{pmatrix} 1 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 1 \end{pmatrix} \begin{pmatrix} U_1 \\ U_2 \\ U_3 \end{pmatrix} = \frac{f_0}{6} \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}$$

The global stiffness matrix K is assembled from all 15 element matrices and the system KU = f is solved with Gaussian elimination with partial pivoting.

### Run

```bash
pip install numpy matplotlib
python poisson-equation/poisson_fem.py
```

---

## 2. Spring-Mass-Damper System — C

### Problem

Find the initial velocity of a damped spring-mass system such that the first full oscillation cycle lasts 6.4 s, given mass m = 1 g, spring constant k = 1 g/s², elongation L = 1 cm:

$$\ddot{x} + 0.8\dot{x} + x = 0, \qquad x(0) = 1 \text{ cm}, \quad \dot{x}(t_{i+1}) = 0$$

The analytical solution is x(t) = e^{-0.4t} cos(0.9165t).

### Rayleigh-Ritz Method

Multiplying by a test function w and integrating over a time element [tᵢ, tᵢ₊₁], then integrating by parts:

$$\int_{t_i}^{t_{i+1}} \left(0.8w\frac{dx}{dt} - \frac{dw}{dt}\frac{dx}{dt} + wx\right)dt = -w\left(\frac{dx}{dt}\right)_{t_i}$$

Using linear shape functions N₁(ξ) = (1-ξ)/2 and N₂(ξ) = (1+ξ)/2 and evaluating with h = 6.4/(n-1) leads to the master element.

### Master Element

$$\begin{pmatrix} 1.5796 & 1.6239 \\ 0.8239 & 2.3796 \end{pmatrix} \begin{pmatrix} u_1 \\ u_2 \end{pmatrix} = \begin{pmatrix} \left(\frac{dx}{dt}\right)_{x_i} \\ \left(\frac{dx}{dt}\right)_{x_{i+1}} = 0 \end{pmatrix}$$

The global system is assembled by overlapping consecutive element matrices and solved with Gauss-Jordan elimination.

### Run

```bash
gcc spring-mass-damper/spring_mass_damper.c -o spring_mass_damper
./spring_mass_damper
```

---

## Repository Structure

```
finite-element-method-/
├── poisson-equation/
│   └── poisson_fem.py        # 2D Poisson equation, linear triangular elements
├── spring-mass-damper/
│   └── spring_mass_damper.c  # Damped oscillator, FEM in time domain
└── README.md
```
