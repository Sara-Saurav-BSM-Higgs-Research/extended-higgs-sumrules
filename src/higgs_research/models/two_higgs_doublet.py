"""Symbolic field definitions for the CP-conserving two-Higgs-doublet model."""

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class HiggsDoublet:
    """A complex SU(2) scalar doublet with real component fields."""

    vev: sp.Symbol
    charged_real: sp.Symbol
    charged_imag: sp.Symbol
    neutral_real: sp.Symbol
    neutral_imag: sp.Symbol

    def as_matrix(self) -> sp.Matrix:
        """Return the complex scalar doublet in canonical normalization."""
        return sp.Matrix(
            [
                (self.charged_real + sp.I * self.charged_imag) / sp.sqrt(2),
                (self.vev + self.neutral_real + sp.I * self.neutral_imag) / sp.sqrt(2),
            ]
        )


def make_higgs_doublet(index: int) -> HiggsDoublet:
    """Construct a symbolic Higgs doublet with real-valued components."""
    if index not in (1, 2):
        raise ValueError("The 2HDM contains doublets indexed by 1 and 2.")

    return HiggsDoublet(
        vev=sp.Symbol(f"v{index}", real=True),
        charged_real=sp.Symbol(f"x{index}", real=True),
        charged_imag=sp.Symbol(f"y{index}", real=True),
        neutral_real=sp.Symbol(f"rho{index}", real=True),
        neutral_imag=sp.Symbol(f"eta{index}", real=True),
    )


def gauge_invariant_bilinear(left: HiggsDoublet, right: HiggsDoublet) -> sp.Expr:
    """Compute Phi_left^dagger Phi_right."""
    return sp.expand((left.as_matrix().conjugate().T * right.as_matrix())[0])


@dataclass(frozen=True)
class TwoHDMParameters:
    """Real parameters of the CP-conserving, softly Z2-broken 2HDM."""

    m11_sq: sp.Symbol
    m22_sq: sp.Symbol
    m12_sq: sp.Symbol
    lambda1: sp.Symbol
    lambda2: sp.Symbol
    lambda3: sp.Symbol
    lambda4: sp.Symbol
    lambda5: sp.Symbol


def make_2hdm_parameters() -> TwoHDMParameters:
    """Construct independent, real-valued potential parameters."""
    names = (
        "m11_sq",
        "m22_sq",
        "m12_sq",
        "lambda1",
        "lambda2",
        "lambda3",
        "lambda4",
        "lambda5",
    )
    symbols = [sp.Symbol(name, real=True) for name in names]
    return TwoHDMParameters(*symbols)


def scalar_potential(
    phi1: HiggsDoublet,
    phi2: HiggsDoublet,
    params: TwoHDMParameters,
) -> sp.Expr:
    """Build the CP-conserving 2HDM potential from field bilinears."""
    b11 = gauge_invariant_bilinear(phi1, phi1)
    b22 = gauge_invariant_bilinear(phi2, phi2)
    b12 = gauge_invariant_bilinear(phi1, phi2)
    b21 = gauge_invariant_bilinear(phi2, phi1)

    return (
        params.m11_sq * b11
        + params.m22_sq * b22
        - params.m12_sq * (b12 + b21)
        + sp.Rational(1, 2) * params.lambda1 * b11**2
        + sp.Rational(1, 2) * params.lambda2 * b22**2
        + params.lambda3 * b11 * b22
        + params.lambda4 * b12 * b21
        + sp.Rational(1, 2) * params.lambda5 * (b12**2 + b21**2)
    )


def neutral_tadpoles(
    phi1: HiggsDoublet,
    phi2: HiggsDoublet,
    params: TwoHDMParameters,
) -> tuple[sp.Expr, sp.Expr]:
    """Return dV/drho1 and dV/drho2 at zero field fluctuations.

    The real VEVs remain arbitrary, including zero; no stationarity conditions
    are imposed and no division by a VEV is performed.
    """
    potential = scalar_potential(phi1, phi2, params)
    vacuum = {
        field: 0
        for phi in (phi1, phi2)
        for field in (
            phi.charged_real,
            phi.charged_imag,
            phi.neutral_real,
            phi.neutral_imag,
        )
    }
    return (
        sp.expand(sp.diff(potential, phi1.neutral_real).subs(vacuum)),
        sp.expand(sp.diff(potential, phi2.neutral_real).subs(vacuum)),
    )
