# Two-Higgs-Doublet Model: Conventions

## Scope

We study the renormalizable, CP-conserving two-Higgs-doublet model (2HDM) with a softly broken $Z_2$ symmetry. The hard-breaking quartic parameters $\lambda_6$ and $\lambda_7$ are set to zero.

This document specifies the conventions used by the symbolic Python implementation.

## Gauge and charge conventions

The gauge group is $SU(2)_L \times U(1)_Y$.

Electric charge is defined by

$$
Q = T_3 + \frac{Y}{2}.
$$

Both Higgs doublets transform as $SU(2)_L$ doublets with hypercharge $Y=1$.

**Alternative notation:** In the convention $Q=T_3+Y$, the same doublets have $Y=1/2$. These notations must not be mixed.

## Field parametrization

The scalar doublets are

$$
\Phi_i =
\begin{pmatrix}
\phi_i^+\\
(v_i+\rho_i+i\eta_i)/\sqrt{2}
\end{pmatrix},
\qquad i=1,2.
$$

For symbolic computations, the charged complex fields are decomposed as

$$
\phi_i^+=(x_i+i y_i)/\sqrt{2},
\qquad
\phi_i^-=(x_i-i y_i)/\sqrt{2}.
$$

The fields $x_i,y_i,\rho_i,\eta_i$ are real.

We take the vacuum expectation values $v_1,v_2$ to be real, with

$$
v^2=v_1^2+v_2^2.
$$

No tadpole equations are imposed at the stage of constructing the potential.

## Gauge-invariant bilinears

Define

$$
B_{ij}=\Phi_i^\dagger\Phi_j.
$$

They satisfy $B_{ji}=B_{ij}^{*}$.

## Scalar potential

The potential is

$$
\begin{aligned}
V={}&m_{11}^2 B_{11}+m_{22}^2 B_{22}
-m_{12}^2(B_{12}+B_{21})\\
&+\frac{\lambda_1}{2}B_{11}^2
+\frac{\lambda_2}{2}B_{22}^2
+\lambda_3 B_{11}B_{22}\\
&+\lambda_4 B_{12}B_{21}
+\frac{\lambda_5}{2}(B_{12}^2+B_{21}^2).
\end{aligned}
$$

All mass parameters and quartic couplings are real.

The Lagrangian contains the potential as $\mathcal{L}\supset -V$.

## Discrete symmetry

The intended $Z_2$ transformation is

$$
\Phi_1\to\Phi_1,\qquad \Phi_2\to-\Phi_2.
$$

All quadratic diagonal terms and all quartic terms are invariant.

The term proportional to $m_{12}^2$ breaks this symmetry softly. The full potential is $Z_2$-invariant when $m_{12}^2=0$.

When testing the transformation in component form, all components of the second doublet—including its vacuum parameter $v_2$—must transform together. Holding $v_2$ fixed describes a different transformation of the shifted fields.

## Validation strategy

The symbolic implementation will be checked using:

1. Exact conjugation identities for the gauge-invariant bilinears.
2. Hermiticity of the full potential.
3. Explicit $Z_2$ transformation of the potential.
4. Special-case limits and independent reconstruction from the field components.

Later tasks will derive stationary-point conditions and mass matrices from the same potential, without inserting expected results into the implementation.

## Primary research references

- Saurav Bania, *Sum Rules in Extended Higgs Sector Models and Their Phenomenological Applications*, Chapter 3.
- Sara Chopra, *Sum Rules for the Extended Higgs Sector Models: A Comparative Analysis of 2HDM, GM and the Septet Model*.

The implementation must independently check equations taken from either thesis and document discrepancies rather than silently modifying conventions.