"""Basic installation and symbolic-computation smoke tests."""

import sympy as sp


def test_symbolic_identity() -> None:
    """Check that SymPy simplifies an exact trigonometric identity."""
    theta = sp.Symbol("theta", real=True)

    expression = sp.sin(theta) ** 2 + sp.cos(theta) ** 2

    assert sp.simplify(expression - 1) == 0


def test_package_import() -> None:
    """Check that the local research package is importable."""
    import higgs_research

    assert higgs_research.__file__ is not None
