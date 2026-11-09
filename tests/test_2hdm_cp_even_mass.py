"""Independent exact checks of the CP-even neutral-vacuum Hessian."""

from fractions import Fraction

import pytest
import sympy as sp

from higgs_research.models.two_higgs_doublet import (
    HiggsDoublet,
    TwoHDMParameters,
    cp_even_mass_matrix,
    make_2hdm_parameters,
    make_higgs_doublet,
    neutral_tadpoles,
    scalar_potential,
)

type Model = tuple[HiggsDoublet, HiggsDoublet, TwoHDMParameters]


@pytest.fixture
def model() -> Model:
    return make_higgs_doublet(1), make_higgs_doublet(2), make_2hdm_parameters()


def component_potential(
    a: HiggsDoublet, b: HiggsDoublet, p: TwoHDMParameters
) -> sp.Expr:
    """Expand real components directly, without production field/bilinear helpers."""
    u, w = a.vev + a.neutral_real, b.vev + b.neutral_real
    n1 = (a.charged_real**2 + a.charged_imag**2 + u**2 + a.neutral_imag**2) / 2
    n2 = (b.charged_real**2 + b.charged_imag**2 + w**2 + b.neutral_imag**2) / 2
    real = (
        a.charged_real * b.charged_real
        + a.charged_imag * b.charged_imag
        + u * w
        + a.neutral_imag * b.neutral_imag
    ) / 2
    imag = (
        a.charged_real * b.charged_imag
        - a.charged_imag * b.charged_real
        + u * b.neutral_imag
        - a.neutral_imag * w
    ) / 2
    return sp.expand(
        p.m11_sq * n1
        + p.m22_sq * n2
        - 2 * p.m12_sq * real
        + p.lambda1 * n1**2 / 2
        + p.lambda2 * n2**2 / 2
        + p.lambda3 * n1 * n2
        + p.lambda4 * (real**2 + imag**2)
        + p.lambda5 * (real**2 - imag**2)
    )


def test_general_coefficients_symmetry_and_immutability(model: Model) -> None:
    a, b, p = model
    v, w = a.vev, b.vev
    matrix = cp_even_mass_matrix(a, b, p)
    assert isinstance(matrix, sp.ImmutableMatrix)
    assert matrix.shape == (2, 2)
    assert matrix == matrix.T
    with pytest.raises(TypeError):
        matrix[0, 0] = 0
    mixed = (p.lambda3 + p.lambda4 + p.lambda5) / 2
    assert sp.Poly(matrix[0, 0], v, w).as_dict() == {
        (0, 0): p.m11_sq,
        (2, 0): 3 * p.lambda1 / 2,
        (0, 2): mixed,
    }
    assert sp.Poly(matrix[1, 1], v, w).as_dict() == {
        (0, 0): p.m22_sq,
        (0, 2): 3 * p.lambda2 / 2,
        (2, 0): mixed,
    }
    assert sp.Poly(matrix[0, 1], v, w).as_dict() == {
        (0, 0): -p.m12_sq,
        (1, 1): 2 * mixed,
    }
    assert matrix.free_symbols <= {v, w, *vars(p).values()}
    assert all(entry.is_polynomial(v, w) for entry in matrix)


def test_independent_component_quadratic_expansion(model: Model) -> None:
    a, b, p = model
    potential = component_potential(a, b, p)
    assert sp.expand(potential - scalar_potential(a, b, p)) == 0
    neutral = potential.subs(
        {
            a.charged_real: 0,
            a.charged_imag: 0,
            a.neutral_imag: 0,
            b.charged_real: 0,
            b.charged_imag: 0,
            b.neutral_imag: 0,
        }
    )
    r, s = a.neutral_real, b.neutral_real
    polynomial = sp.Poly(neutral, r, s)
    matrix = cp_even_mass_matrix(a, b, p)
    assert sp.expand(2 * polynomial.coeff_monomial(r**2) - matrix[0, 0]) == 0
    assert sp.expand(2 * polynomial.coeff_monomial(s**2) - matrix[1, 1]) == 0
    assert sp.expand(polynomial.coeff_monomial(r * s) - matrix[0, 1]) == 0
    quadratic = sum(c * r**i * s**j for (i, j), c in polynomial.terms() if i + j == 2)
    rho = sp.Matrix([r, s])
    assert sp.expand(quadratic - (rho.T * matrix * rho)[0] / 2) == 0
    v0 = neutral.subs({r: 0, s: 0})
    assert (sp.hessian(v0, (a.vev, b.vev)) - matrix).applyfunc(sp.expand).is_zero_matrix


def test_stationarity_recipe_and_division_free_identities(model: Model) -> None:
    a, b, p = model
    v, w = a.vev, b.vev
    matrix = cp_even_mass_matrix(a, b, p)
    t1, t2 = neutral_tadpoles(a, b, p)
    assert sp.expand(v * (matrix[0, 0] - p.lambda1 * v**2) - p.m12_sq * w - t1) == 0
    assert sp.expand(w * (matrix[1, 1] - p.lambda2 * w**2) - p.m12_sq * v - t2) == 0
    # Only this recipe assumes BOTH VEVs nonzero, before solving/dividing.
    vn, wn = sp.symbols("vn wn", real=True, nonzero=True)
    vevs = {v: vn, w: wn}
    tadpoles = [t.subs(vevs) for t in (t1, t2)]
    solution = sp.solve(tadpoles, (p.m11_sq, p.m22_sq), dict=True)[0]
    assert all(sp.simplify(t.subs(solution, simultaneous=True)) == 0 for t in tadpoles)
    stationary = matrix.subs(vevs).subs(solution, simultaneous=True)
    expected = sp.Matrix(
        [
            [
                p.lambda1 * vn**2 + p.m12_sq * wn / vn,
                -p.m12_sq + (p.lambda3 + p.lambda4 + p.lambda5) * vn * wn,
            ],
            [
                -p.m12_sq + (p.lambda3 + p.lambda4 + p.lambda5) * vn * wn,
                p.lambda2 * wn**2 + p.m12_sq * vn / wn,
            ],
        ]
    )
    assert (stationary - expected).applyfunc(sp.simplify).is_zero_matrix


@pytest.mark.parametrize("zero_index", [0, 1])
def test_single_zero_vev_stationary_branches(model: Model, zero_index: int) -> None:
    a, b, p = model
    vevs = (a.vev, b.vev)
    masses = (p.m11_sq, p.m22_sq)
    quartics = (p.lambda1, p.lambda2)
    active = 1 - zero_index
    nonzero = sp.Symbol("u", real=True, nonzero=True)
    branch = {vevs[zero_index]: 0, vevs[active]: nonzero}
    tadpoles = [t.subs(branch) for t in neutral_tadpoles(a, b, p)]
    solution = sp.solve(tadpoles, (p.m12_sq, masses[active]), dict=True)[0]
    assert solution[p.m12_sq] == 0
    assert solution[masses[active]] == -quartics[active] * nonzero**2 / 2
    matrix = cp_even_mass_matrix(a, b, p).subs(branch)
    # The general, nonstationary branch retains soft mixing.
    assert matrix[0, 1] == -p.m12_sq
    stationary = matrix.subs(solution, simultaneous=True)
    assert stationary[0, 1] == 0
    assert stationary[active, active] == quartics[active] * nonzero**2
    assert (
        sp.expand(
            stationary[zero_index, zero_index]
            - masses[zero_index]
            - (p.lambda3 + p.lambda4 + p.lambda5) * nonzero**2 / 2
        )
        == 0
    )
    assert all(sp.simplify(t.subs(solution, simultaneous=True)) == 0 for t in tadpoles)


def test_origin_and_special_limits(model: Model) -> None:
    a, b, p = model
    matrix = cp_even_mass_matrix(a, b, p)
    quadratic = sp.Matrix([[p.m11_sq, -p.m12_sq], [-p.m12_sq, p.m22_sq]])
    origin = {a.vev: 0, b.vev: 0}
    assert matrix.subs(origin) == quadratic
    assert all(t.subs(origin) == 0 for t in neutral_tadpoles(a, b, p))
    quartics = dict.fromkeys((p.lambda1, p.lambda2, p.lambda3, p.lambda4, p.lambda5), 0)
    assert matrix.subs(quartics) == quadratic
    assert matrix.subs(quartics | {p.m11_sq: 0, p.m22_sq: 0}) == sp.Matrix(
        [[0, -p.m12_sq], [-p.m12_sq, 0]]
    )
    decoupled = {p.m12_sq: 0, p.lambda3: 0, p.lambda4: 0, p.lambda5: 0}
    assert matrix.subs(decoupled) == sp.diag(
        p.m11_sq + 3 * p.lambda1 * a.vev**2 / 2, p.m22_sq + 3 * p.lambda2 * b.vev**2 / 2
    )
    cancel_mixing = {p.m12_sq: (p.lambda3 + p.lambda4 + p.lambda5) * a.vev * b.vev}
    assert sp.expand(matrix.subs(cancel_mixing)[0, 1]) == 0


def test_symbolic_invariants_and_eigenvalues(model: Model) -> None:
    a, b, p = model
    matrix = cp_even_mass_matrix(a, b, p)
    # Independent coefficient extraction supplies the invariant oracle.
    neutral = component_potential(a, b, p).subs(
        {
            a.charged_real: 0,
            a.charged_imag: 0,
            a.neutral_imag: 0,
            b.charged_real: 0,
            b.charged_imag: 0,
            b.neutral_imag: 0,
        }
    )
    poly = sp.Poly(neutral, a.neutral_real, b.neutral_real)
    aa = 2 * poly.coeff_monomial(a.neutral_real**2)
    dd = 2 * poly.coeff_monomial(b.neutral_real**2)
    bb = poly.coeff_monomial(a.neutral_real * b.neutral_real)
    trace, determinant = aa + dd, aa * dd - bb**2
    assert sp.expand(matrix.trace() - trace) == 0
    assert sp.expand(matrix.det() - determinant) == 0
    z = sp.Symbol("z")
    assert (
        sp.expand(matrix.charpoly(z).as_expr() - (z**2 - trace * z + determinant)) == 0
    )
    discriminant = sp.sqrt((aa - dd) ** 2 + 4 * bb**2)
    low, high = (trace - discriminant) / 2, (trace + discriminant) / 2
    assert sp.expand(low + high - trace) == 0
    assert sp.expand(low * high - determinant) == 0
    for eigenvalue in (low, high):
        assert (
            sp.simplify(z**2 - trace * z + determinant).subs(z, eigenvalue).expand()
            == 0
        )


def test_distinct_and_degenerate_orthonormal_eigenvectors(model: Model) -> None:
    a, b, p = model
    # Stationary examples: equal VEVs, lambda345=1, lambda1=lambda2=2.
    common = {
        a.vev: 1,
        b.vev: 1,
        p.lambda1: 2,
        p.lambda2: 2,
        p.lambda3: 1,
        p.lambda4: 0,
        p.lambda5: 0,
    }
    distinct = common | {p.m12_sq: sp.Rational(1, 2), p.m11_sq: -1, p.m22_sq: -1}
    degenerate = common | {
        p.m12_sq: 1,
        p.m11_sq: -sp.Rational(1, 2),
        p.m22_sq: -sp.Rational(1, 2),
    }
    rotation = sp.Matrix([[1, 1], [-1, 1]]) / sp.sqrt(2)
    assert rotation.T * rotation == sp.eye(2)
    for values, expected in ((distinct, sp.diag(2, 3)), (degenerate, 3 * sp.eye(2))):
        assert all(t.subs(values) == 0 for t in neutral_tadpoles(a, b, p))
        matrix = cp_even_mass_matrix(a, b, p).subs(values)
        assert matrix.eigenvals() == (
            {sp.Integer(2): 1, sp.Integer(3): 1}
            if values is distinct
            else {sp.Integer(3): 2}
        )
        assert rotation.T * matrix * rotation == expected
        for i in range(2):
            assert matrix * rotation[:, i] == expected[i, i] * rotation[:, i]


def test_fixed_exact_rational_benchmark(model: Model) -> None:
    a, b, p = model
    values = {
        a.vev: 2,
        b.vev: 3,
        p.m11_sq: 5,
        p.m22_sq: 7,
        p.m12_sq: 11,
        p.lambda1: sp.Rational(1, 2),
        p.lambda2: sp.Rational(3, 4),
        p.lambda3: sp.Rational(1, 5),
        p.lambda4: -sp.Rational(2, 5),
        p.lambda5: sp.Rational(3, 5),
    }
    matrix = cp_even_mass_matrix(a, b, p).subs(values)
    expected = sp.ImmutableMatrix(
        [
            [sp.Rational(49, 5), -sp.Rational(43, 5)],
            [-sp.Rational(43, 5), sp.Rational(717, 40)],
        ]
    )
    assert matrix == expected
    assert matrix.trace() == sp.Rational(1109, 40)
    assert matrix.det() == sp.Rational(20341, 200)

    # Independent Fraction-only potential evaluations and exact finite differences.
    # Richardson cancellation removes the quartic diagonal finite-step error;
    # the mixed central difference is exact for this biquadratic potential.
    def potential(x: int, y: int) -> Fraction:
        u, w = Fraction(2 + x), Fraction(3 + y)
        return (
            5 * u * u / 2
            + 7 * w * w / 2
            - 11 * u * w
            + u**4 / 16
            + 3 * w**4 / 32
            + u * u * w * w / 10
        )

    def diagonal(axis: int, step: int) -> Fraction:
        plus = (step, 0) if axis == 0 else (0, step)
        minus = (-step, 0) if axis == 0 else (0, -step)
        return (potential(*plus) - 2 * potential(0, 0) + potential(*minus)) / step**2

    assert (4 * diagonal(0, 1) - diagonal(0, 2)) / 3 == Fraction(49, 5)
    assert (4 * diagonal(1, 1) - diagonal(1, 2)) / 3 == Fraction(717, 40)
    mixed = (
        potential(1, 1) - potential(1, -1) - potential(-1, 1) + potential(-1, -1)
    ) / 4
    assert mixed == Fraction(-43, 5)
