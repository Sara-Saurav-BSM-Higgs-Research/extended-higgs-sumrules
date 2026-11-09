"""Tests for the symbolic 2HDM field representation."""

import sympy as sp

from higgs_research.models.two_higgs_doublet import (
    gauge_invariant_bilinear,
    make_2hdm_parameters,
    make_higgs_doublet,
    scalar_potential,
)


def test_doublet_has_two_components() -> None:
    """Each SU(2) scalar doublet must contain two complex components."""
    doublet = make_higgs_doublet(1)

    assert doublet.as_matrix().shape == (2, 1)


def test_bilinear_conjugation() -> None:
    """Verify B_21 = conjugate(B_12) directly from field components."""
    phi1 = make_higgs_doublet(1)
    phi2 = make_higgs_doublet(2)

    b12 = gauge_invariant_bilinear(phi1, phi2)
    b21 = gauge_invariant_bilinear(phi2, phi1)

    assert sp.simplify(b21 - sp.conjugate(b12)) == 0


def test_diagonal_bilinear_is_real() -> None:
    """Verify that Phi_i^dagger Phi_i is real."""
    phi1 = make_higgs_doublet(1)
    b11 = gauge_invariant_bilinear(phi1, phi1)

    assert sp.simplify(b11 - sp.conjugate(b11)) == 0


def test_invalid_doublet_index() -> None:
    """Only doublets 1 and 2 are valid in this model."""
    import pytest

    with pytest.raises(ValueError):
        make_higgs_doublet(3)


def test_scalar_potential_is_real() -> None:
    """The CP-conserving potential must be Hermitian."""
    phi1 = make_higgs_doublet(1)
    phi2 = make_higgs_doublet(2)
    params = make_2hdm_parameters()

    potential = scalar_potential(phi1, phi2, params)

    assert sp.expand(potential - sp.conjugate(potential)) == 0


def test_quartic_free_limit() -> None:
    """Switching off quartics must leave the quadratic potential."""
    phi1 = make_higgs_doublet(1)
    phi2 = make_higgs_doublet(2)
    params = make_2hdm_parameters()

    potential = scalar_potential(phi1, phi2, params)

    zero_quartics = {
        params.lambda1: 0,
        params.lambda2: 0,
        params.lambda3: 0,
        params.lambda4: 0,
        params.lambda5: 0,
    }

    expected = (
        params.m11_sq * gauge_invariant_bilinear(phi1, phi1)
        + params.m22_sq * gauge_invariant_bilinear(phi2, phi2)
        - params.m12_sq
        * (gauge_invariant_bilinear(phi1, phi2) + gauge_invariant_bilinear(phi2, phi1))
    )

    assert sp.expand(potential.subs(zero_quartics) - expected) == 0


def test_z2_soft_breaking_identity() -> None:
    """Only the soft m12_sq term changes under Phi2 -> -Phi2."""
    phi1 = make_higgs_doublet(1)
    phi2 = make_higgs_doublet(2)
    params = make_2hdm_parameters()

    potential = scalar_potential(phi1, phi2, params)

    # Flip the entire second doublet, including its vacuum parameter.
    z2_substitution = {
        phi2.vev: -phi2.vev,
        phi2.charged_real: -phi2.charged_real,
        phi2.charged_imag: -phi2.charged_imag,
        phi2.neutral_real: -phi2.neutral_real,
        phi2.neutral_imag: -phi2.neutral_imag,
    }

    transformed = potential.xreplace(z2_substitution)

    b12 = gauge_invariant_bilinear(phi1, phi2)
    b21 = gauge_invariant_bilinear(phi2, phi1)

    expected_difference = 2 * params.m12_sq * (b12 + b21)

    assert sp.expand(transformed - potential - expected_difference) == 0


def test_exact_z2_symmetry_without_soft_breaking() -> None:
    """The complete potential is Z2 invariant when m12_sq vanishes."""
    phi1 = make_higgs_doublet(1)
    phi2 = make_higgs_doublet(2)
    params = make_2hdm_parameters()

    potential = scalar_potential(phi1, phi2, params)
    potential_without_soft_breaking = potential.subs(params.m12_sq, 0)

    z2_substitution = {
        phi2.vev: -phi2.vev,
        phi2.charged_real: -phi2.charged_real,
        phi2.charged_imag: -phi2.charged_imag,
        phi2.neutral_real: -phi2.neutral_real,
        phi2.neutral_imag: -phi2.neutral_imag,
    }

    transformed = potential_without_soft_breaking.xreplace(z2_substitution)

    assert sp.expand(transformed - potential_without_soft_breaking) == 0


def test_neutral_vacuum_potential() -> None:
    """Check the potential evaluated at the neutral field configuration."""
    phi1 = make_higgs_doublet(1)
    phi2 = make_higgs_doublet(2)
    params = make_2hdm_parameters()

    potential = scalar_potential(phi1, phi2, params)

    vacuum_substitution = {
        phi1.charged_real: 0,
        phi1.charged_imag: 0,
        phi1.neutral_real: 0,
        phi1.neutral_imag: 0,
        phi2.charged_real: 0,
        phi2.charged_imag: 0,
        phi2.neutral_real: 0,
        phi2.neutral_imag: 0,
    }

    vacuum_potential = potential.subs(vacuum_substitution)

    lambda345 = params.lambda3 + params.lambda4 + params.lambda5

    expected = (
        params.m11_sq * phi1.vev**2 / 2
        + params.m22_sq * phi2.vev**2 / 2
        - params.m12_sq * phi1.vev * phi2.vev
        + params.lambda1 * phi1.vev**4 / 8
        + params.lambda2 * phi2.vev**4 / 8
        + lambda345 * phi1.vev**2 * phi2.vev**2 / 4
    )

    assert sp.expand(vacuum_potential - expected) == 0


def test_imaginary_bilinear_distinguishes_lambda4_and_lambda5() -> None:
    """For Phi1=(1,0)/sqrt(2), Phi2=(i,0)/sqrt(2), B12=i/2."""
    phi1 = make_higgs_doublet(1)
    phi2 = make_higgs_doublet(2)
    params = make_2hdm_parameters()
    configuration = {
        phi1.vev: 0,
        phi1.charged_real: 1,
        phi1.charged_imag: 0,
        phi1.neutral_real: 0,
        phi1.neutral_imag: 0,
        phi2.vev: 0,
        phi2.charged_real: 0,
        phi2.charged_imag: 1,
        phi2.neutral_real: 0,
        phi2.neutral_imag: 0,
        params.m11_sq: 0,
        params.m22_sq: 0,
        params.m12_sq: 0,
        params.lambda1: 0,
        params.lambda2: 0,
        params.lambda3: 0,
    }
    assert gauge_invariant_bilinear(phi1, phi2).subs(configuration) == sp.I / 2

    actual = sp.expand(scalar_potential(phi1, phi2, params).subs(configuration))
    # B12*B21=1/4, whereas (B12**2+B21**2)/2=-1/4.
    expected = (params.lambda4 - params.lambda5) / 4
    assert sp.expand(actual - expected) == 0
    exchanged = actual.xreplace(
        {params.lambda4: params.lambda5, params.lambda5: params.lambda4}
    )
    assert sp.expand(exchanged - expected) != 0
