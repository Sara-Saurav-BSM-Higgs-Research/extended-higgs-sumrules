# Extended Higgs Sector Sum Rules and Phenomenology

**An open-source computational framework for symbolic derivations, perturbative unitarity constraints, and phenomenological studies of extended Higgs sectors.**

**Authors:** Saurav Bania and Sara Chopra

**Research organization:** [Sara-Saurav-BSM-Higgs-Research](https://github.com/Sara-Saurav-BSM-Higgs-Research)

**Project status:** Initial development — infrastructure setup and analytical reconstruction.

---

## Overview

The scalar sector of the Standard Model (SM) provides the mechanism of electroweak symmetry breaking and generates the masses of the electroweak gauge bosons. Extensions of this sector introduce additional scalar degrees of freedom and modified Higgs–gauge-boson interactions, subject to theoretical consistency conditions and experimental constraints.

One important consistency requirement arises from perturbative unitarity in longitudinal vector-boson scattering. At high center-of-mass energies, contributions from gauge and scalar interactions must combine to cancel terms exhibiting unphysical energy growth. These cancellations lead to relations among scalar–gauge-boson couplings known as **unitarity sum rules**.

This repository aims to develop a reproducible computational framework for deriving, verifying and studying these relations in extended Higgs sector models.

The project builds upon two master's theses completed at the University of Hamburg:

- **Saurav Bania:** *Sum Rules in Extended Higgs Sector Models and Their Phenomenological Applications* — analytical studies of the Two-Higgs-Doublet Model (2HDM) and Georgi–Machacek (GM) model, including gauge-boson scattering, coupling sum rules, perturbative unitarity bounds and fermionic interactions.

- **Sara Chopra:** *Sum Rules for the Extended Higgs Sector Models: A Comparative Analysis of 2HDM, GM and the Septet Model* — comparative analytical and phenomenological investigations of the 2HDM, GM and scalar Septet models, incorporating theoretical constraints, Higgs coupling limits and numerical studies using HiggsTools.

The initial objective is to reconstruct and independently validate the analytical calculations through a documented, testable, Python-based implementation. Subsequent development will incorporate and extend the numerical analysis, including studies with updated experimental datasets.

## Scientific Objectives

The project is organized around five principal objectives.

### 1. Symbolic reconstruction of extended Higgs sectors

Implement the scalar field content, potentials and electroweak interactions of the models under investigation.

The intended calculations include:

- Scalar potentials and vacuum configurations
- Tadpole and minimization conditions
- CP-even, CP-odd and charged scalar mass matrices
- Mass eigenstates, mixing transformations and Goldstone modes
- Higgs interactions with electroweak gauge bosons and fermions

### 2. High-energy scattering amplitudes and sum rules

Reconstruct tree-level scattering amplitudes involving longitudinally polarized electroweak gauge bosons, including

$$
W_L^+W_L^- \rightarrow W_L^+W_L^-,
$$

$$
W_L^+W_L^- \rightarrow Z_LZ_L,
$$

$$
Z_LZ_L \rightarrow Z_LZ_L.
$$

The symbolic workflow will construct relevant contributions, perform Lorentz contractions, implement scattering kinematics, extract high-energy behavior and derive the conditions required for cancellation of energy-growing terms.

The 2HDM provides a useful reference identity:

$$
(\kappa_V^h)^2+(\kappa_V^H)^2=1,
$$

where the normalized gauge-boson couplings are

$$
\kappa_V^h=\sin(\beta-\alpha),
\qquad
\kappa_V^H=\cos(\beta-\alpha).
$$

Models containing higher electroweak representations have additional scalar contributions to vector-boson scattering. Their sum rules will be reconstructed using the appropriate model-specific coupling conventions and scalar spectra.

The objective is to **derive these relations from the underlying interactions**, rather than merely verify identities supplied as inputs.

### 3. Perturbative unitarity and theoretical constraints

Investigate the implications of partial-wave unitarity for scalar masses, mixing angles and interaction strengths.

For a conventional elastic scattering partial-wave normalization, the tree-level perturbative criterion is

$$
\left|\operatorname{Re} a_0\right| \leq \frac{1}{2}.
$$

Planned investigations include individual scattering channels, mass and coupling bounds, and extensions to coupled-channel analyses.

### 4. Comparison of extended scalar sectors

Develop a common computational framework for:

| Model | Scalar extension | Key focus |
|---|---|---|
| Two-Higgs-Doublet Model (2HDM) | Additional complex electroweak doublet | Scalar mixing, gauge couplings and unitarity sum rules |
| Georgi–Machacek model (GM) | Real and complex electroweak triplets | Custodial multiplets, fiveplet interactions and unitarity constraints |
| Scalar Septet model | Higher-isospin electroweak scalar representation | Extended charge spectrum and modified scattering sum rules |

The initial implementation will focus on the 2HDM, followed by the GM model and the Septet extension.

### 5. Phenomenology and collider constraints

Integrate analytical model predictions with numerical parameter studies and experimental constraints.

Future development is intended to include:

- Reproduction of existing numerical results
- Integration of the original Python analysis
- HiggsTools-based exclusion studies
- Comparison of theoretical and experimentally allowed parameter regions
- Evaluation of updated HiggsBounds datasets and applicable mass ranges
- Investigation of additional phenomenological scenarios

Experimental datasets, software releases and model assumptions will be versioned and documented to support meaningful comparisons.

## Computational Approach

The primary implementation will use **Python 3 and SymPy** for symbolic calculations, supported by numerical scientific Python tools.

The software strategy is designed to avoid mandatory dependence on proprietary computational systems.

| Component | Intended technology |
|---|---|
| Symbolic mathematics | Python, SymPy |
| Numerical computation | NumPy, SciPy |
| Scientific visualization | Matplotlib |
| Interactive derivations | Jupyter notebooks, VS Code |
| Automated testing | pytest |
| Code formatting and linting | Ruff |
| Collider phenomenology | HiggsTools |
| Advanced tensor algebra | Cadabra, when required |
| Optional independent verification | Wolfram Language / Mathematica |
| Performance-oriented extensions | C++ or FORM, when justified |

The analytical calculations will be implemented in reusable Python modules. Jupyter notebooks will primarily serve as documented research workflows, demonstrations and exploratory studies.

## Verification and Validation

Scientific reproducibility is a central design requirement.

Three complementary levels of verification and validation are planned.

**Algebraic verification:** Automated checks of symbolic identities, minimization conditions, matrix diagonalization, Goldstone eigenvalues, coupling relations and high-energy cancellations.

**Independent physics validation:** Comparison against independently derived results, published analytical expressions, consistency limits, symmetry relations and, where practical, alternative symbolic implementations.

**Regression testing:** Reproduction of reference calculations, numerical benchmarks and selected published figures, with automated checks to detect unintended changes after software modifications.

All calculations will explicitly document assumptions, field definitions, conventions and normalizations.

Where reconstruction reveals discrepancies between independently derived results and the source theses, these will be investigated and documented rather than silently modified.

## Planned Repository Organization

```text
Masters_Research_SS/
├── src/higgs_research/
│   ├── core/
│   ├── models/
│   ├── scattering/
│   ├── unitarity/
│   └── phenomenology/
├── notebooks/
├── tests/
├── docs/
├── studies/
├── mathematica/
├── data/
├── scripts/
├── pyproject.toml
├── README.md
├── CONTRIBUTING.md
├── CITATION.cff
└── LICENSE
```

This is the intended project layout and will be populated progressively as the implementation develops.

## Development Roadmap

| Phase | Scope | Status |
|---|---|---|
| 0 | Repository, environment, package architecture and testing infrastructure | In progress |
| 1 | 2HDM potential, tadpoles, mass matrices and physical states | Planned |
| 2 | 2HDM electroweak interactions and Higgs couplings | Planned |
| 3 | Longitudinal gauge-boson scattering and sum-rule derivations | Planned |
| 4 | Perturbative unitarity and theoretical mass bounds | Planned |
| 5 | GM model reconstruction and custodial multiplet analysis | Planned |
| 6 | Septet model and comparative calculations | Planned |
| 7 | Existing Python phenomenology and HiggsTools integration | Planned |
| 8 | Updated experimental analyses and research extensions | Planned |

The roadmap may evolve as the underlying physics and software requirements are evaluated.

## Reproducibility and Open Science

This project aims to follow FAIR research principles: Findability, Accessibility, Interoperability and Reusability.

The development workflow will emphasize:

- Publicly accessible source code
- Explicit software dependencies and environment documentation
- Version-controlled analytical and numerical workflows
- Automated physics verification and software tests
- Transparent scientific assumptions and calculation conventions
- Traceability between analytical derivations and reference equations
- Documented provenance of experimental datasets
- Clear attribution of external software and scientific references

Where third-party datasets or software are used, their respective licences and redistribution conditions will be respected.

Reproducible releases and archival research outputs may be published as the framework matures.

## Getting Started

The package and installation instructions are under development.

Once the initial Python infrastructure is established, this section will contain verified installation instructions, commands for running the automated tests, and examples of the symbolic derivations.

At present, the repository is in the setup phase; the planned scientific capabilities described above should not yet be interpreted as implemented functionality.

## Research Origins and Acknowledgments

This project originates from the master's research of Saurav Bania and Sara Chopra at the University of Hamburg, under the supervision of **Prof. Dr. Gudrid Moortgat-Pick** and **Prof. Dr. Sven Heinemeyer**.

The authors acknowledge the scientific guidance and contributions described in their respective master's theses.

The present repository is an independently developed, collaborative continuation of that research. It is not presented as an official software project of the University of Hamburg, CERN or any other institution.

## Authors and Contributions

**Saurav Bania** — original analytical studies of unitarity sum rules and phenomenological applications in the 2HDM and GM model; collaborative development of the computational framework.

**Sara Chopra** — comparative analytical and phenomenological studies of the 2HDM, GM and Septet models, including HiggsTools-based numerical analyses; collaborative development and extension of the research framework.

Specific code contributions and scientific results will be documented through the project's Git history, pull requests and subsequent research outputs.

## Licence and Citation

This project's original software and documentation are released under the **BSD 3-Clause License**. See the [LICENSE](LICENSE) file for the full terms and conditions.

The licence permits use, modification and redistribution, including commercial use, subject to its conditions. The names of the copyright holders and contributors may not be used to endorse derivative products without prior written permission.

External software, datasets and third-party research materials remain subject to their respective licences and attribution requirements.

Scientific users are encouraged to acknowledge the contributions of **Saurav Bania and Sara Chopra** when using this framework in research.

Machine-readable citation metadata (`CITATION.cff`) and DOI-backed release information will be provided as the project matures.

---

**Project principle:** Reconstruct the physics, verify the mathematics, document the assumptions and make the research reproducible.
