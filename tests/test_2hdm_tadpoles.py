"""Independent neutral-vacuum and exact polynomial checks of 2HDM tadpoles."""

import pytest
import sympy as sp

from higgs_research.models.two_higgs_doublet import (
    HiggsDoublet,
    TwoHDMParameters,
    make_2hdm_parameters,
    make_higgs_doublet,
    neutral_tadpoles,
    scalar_potential,
)


def independent_vacuum_potential(
    v1: sp.Symbol, v2: sp.Symbol, p: TwoHDMParameters
) -> sp.Expr:
    """Reconstruct V0 from Phi_i=(0,v_i/sqrt(2)), without production helpers."""
    return (
        p.m11_sq * v1**2 / 2
        + p.m22_sq * v2**2 / 2
        - p.m12_sq * v1 * v2
        + p.lambda1 * v1**4 / 8
        + p.lambda2 * v2**4 / 8
        + p.lambda3 * v1**2 * v2**2 / 4
        + p.lambda4 * v1**2 * v2**2 / 4
        + p.lambda5 * v1**2 * v2**2 / 4
    )


def zero_fluctuations(*doublets: HiggsDoublet) -> dict[sp.Symbol, int]:
    return {
        field: 0
        for phi in doublets
        for field in (
            phi.charged_real,
            phi.charged_imag,
            phi.neutral_real,
            phi.neutral_imag,
        )
    }


def test_tadpoles_match_independent_vacuum_derivatives() -> None:
    phi1, phi2 = make_higgs_doublet(1), make_higgs_doublet(2)
    p = make_2hdm_parameters()
    v0 = independent_vacuum_potential(phi1.vev, phi2.vev, p)
    actual_v0 = scalar_potential(phi1, phi2, p).subs(zero_fluctuations(phi1, phi2))
    assert sp.expand(actual_v0 - v0) == 0
    tadpoles = neutral_tadpoles(phi1, phi2, p)
    assert isinstance(tadpoles, tuple) and len(tadpoles) == 2
    for tadpole, vev in zip(tadpoles, (phi1.vev, phi2.vev), strict=True):
        assert sp.expand(tadpole - sp.diff(v0, vev)) == 0


def test_exact_tadpole_polynomial_coefficients() -> None:
    phi1, phi2 = make_higgs_doublet(1), make_higgs_doublet(2)
    p = make_2hdm_parameters()
    v1, v2 = phi1.vev, phi2.vev
    t1, t2 = neutral_tadpoles(phi1, phi2, p)
    # Check the complete VEV monomial maps, including absence of extra terms.
    assert sp.Poly(t1, v1, v2).as_dict() == {
        (1, 0): p.m11_sq,
        (0, 1): -p.m12_sq,
        (3, 0): p.lambda1 / 2,
        (1, 2): p.lambda3 / 2 + p.lambda4 / 2 + p.lambda5 / 2,
    }
    assert sp.Poly(t2, v1, v2).as_dict() == {
        (0, 1): p.m22_sq,
        (1, 0): -p.m12_sq,
        (0, 3): p.lambda2 / 2,
        (2, 1): p.lambda3 / 2 + p.lambda4 / 2 + p.lambda5 / 2,
    }
    for coupling in (p.lambda3, p.lambda4, p.lambda5):
        assert t1.coeff(coupling) == v1 * v2**2 / 2
        assert t2.coeff(coupling) == v1**2 * v2 / 2
    assert t1.free_symbols | t2.free_symbols <= {
        v1,
        v2,
        p.m11_sq,
        p.m22_sq,
        p.m12_sq,
        p.lambda1,
        p.lambda2,
        p.lambda3,
        p.lambda4,
        p.lambda5,
    }


@pytest.mark.parametrize("v1,v2", [(0, 0), (0, 3), (2, 0)])
def test_zero_vev_cases(v1: int, v2: int) -> None:
    phi1, phi2 = make_higgs_doublet(1), make_higgs_doublet(2)
    p = make_2hdm_parameters()
    t1, t2 = neutral_tadpoles(phi1, phi2, p)
    vacuum = {phi1.vev: v1, phi2.vev: v2}
    # Mixed quartics vanish when either VEV is zero; the soft term need not.
    assert t1.subs(vacuum) == p.m11_sq * v1 - p.m12_sq * v2 + p.lambda1 * v1**3 / 2
    assert t2.subs(vacuum) == p.m22_sq * v2 - p.m12_sq * v1 + p.lambda2 * v2**3 / 2


def test_quartic_free_and_soft_only_limits() -> None:
    phi1, phi2 = make_higgs_doublet(1), make_higgs_doublet(2)
    p = make_2hdm_parameters()
    t1, t2 = neutral_tadpoles(phi1, phi2, p)
    no_quartics = dict.fromkeys(
        (p.lambda1, p.lambda2, p.lambda3, p.lambda4, p.lambda5), 0
    )
    assert t1.subs(no_quartics) == p.m11_sq * phi1.vev - p.m12_sq * phi2.vev
    assert t2.subs(no_quartics) == p.m22_sq * phi2.vev - p.m12_sq * phi1.vev
    soft_only = no_quartics | {p.m11_sq: 0, p.m22_sq: 0}
    assert t1.subs(soft_only) == -p.m12_sq * phi2.vev
    assert t2.subs(soft_only) == -p.m12_sq * phi1.vev


def test_decoupled_doublet_limit() -> None:
    phi1, phi2 = make_higgs_doublet(1), make_higgs_doublet(2)
    p = make_2hdm_parameters()
    t1, t2 = neutral_tadpoles(phi1, phi2, p)
    decoupled = {p.m12_sq: 0, p.lambda3: 0, p.lambda4: 0, p.lambda5: 0}
    assert t1.subs(decoupled) == p.m11_sq * phi1.vev + p.lambda1 * phi1.vev**3 / 2
    assert t2.subs(decoupled) == p.m22_sq * phi2.vev + p.lambda2 * phi2.vev**3 / 2


def test_exact_rational_benchmark() -> None:
    phi1, phi2 = make_higgs_doublet(1), make_higgs_doublet(2)
    p = make_2hdm_parameters()
    benchmark = {
        phi1.vev: 2,
        phi2.vev: 3,
        p.m11_sq: 5,
        p.m22_sq: 7,
        p.m12_sq: 11,
        p.lambda1: sp.Rational(1, 2),
        p.lambda2: sp.Rational(3, 4),
        p.lambda3: sp.Rational(1, 5),
        p.lambda4: -sp.Rational(2, 5),
        p.lambda5: sp.Rational(3, 5),
    }
    v0 = independent_vacuum_potential(phi1.vev, phi2.vev, p)
    assert v0.subs(benchmark) == -sp.Rational(1969, 160)
    actual_v0 = scalar_potential(phi1, phi2, p).subs(zero_fluctuations(phi1, phi2))
    assert actual_v0.subs(benchmark) == -sp.Rational(1969, 160)
    t1, t2 = neutral_tadpoles(phi1, phi2, p)
    assert t1.subs(benchmark) == -sp.Rational(87, 5)
    assert t2.subs(benchmark) == sp.Rational(461, 40)
