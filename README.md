# MultiContrasts <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)](https://github.com/zhiyuemo/MultiContrasts)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

## Doubly Robust Estimation of Multiple Causal Contrasts

**Author:** [Zhiyue Mo](https://github.com/zhiyuemo)

## What's MultiContrasts?

`MultiContrasts` provides a unified interface for estimating a broad class of causal contrasts using cross-fitted Augmented Inverse Probability Weighting (AIPW) with SuperLearner for nuisance parameter estimation. The package implements the doubly robust, semiparametric efficient estimators described in our accompanying paper.

Rather than restricting investigators to a single effect measure, `MultiContrasts` simultaneously estimates **seven** population-level causal contrasts from a single function call:

| Contrast | Formula | CI Method |
|---|---|---|
| **Risk Difference (RD)** | ψ₁ − ψ₀ | Wald |
| **Risk Ratio (RR)** | ψ₁ / ψ₀ | Log-transform |
| **Odds Ratio (OR)** | [ψ₁/(1−ψ₁)] / [ψ₀/(1−ψ₀)] | Log-transform |
| **Vaccine Efficacy (VE)** | 1 − RR | Via log RR |
| **Excess Relative Risk (ERR)** | (ψ₁ − ψ₀) / ψ₀ | Via log RR |
| **Number Needed to Treat (NNT)** | 1 / \|RD\| | Log-transform or RD flip |
| **Survival Ratio (SR)** | (1−ψ₁) / (1−ψ₀) | Log-transform |

where ψ₁ = E[Y(1)] and ψ₀ = E[Y(0)] are the counterfactual outcome means.

### Why MultiContrasts?

- **Science-first workflow:** Choose the estimand that matches your scientific question, not the one that's easiest to compute.
- **Doubly robust:** Consistent if either the propensity score or outcome model is correctly specified.
- **Rate doubly robust:** Valid root-n inference when nuisance estimators converge at sufficiently fast rates via cross-fitting.
- **Flexible ML:** Separate SuperLearner libraries for propensity score and outcome models.
- **CONSORT 2025 compliant:** Report both absolute and relative measures for binary outcomes from a single analysis.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("zhiyuemo/MultiContrasts")
```

## Quick Start

```r
library(MultiContrasts)
library(survival)

# Prepare data
data("cancer", package = "survival")
cancer <- cancer[!is.na(cancer$ph.ecog) & !is.na(cancer$age) & !is.na(cancer$sex), ]
cancer$dead <- as.integer(cancer$status == 2)
cancer$ecog_bin <- as.integer(cancer$ph.ecog >= 1)
cancer$dead_365 <- ifelse(cancer$time < 365 & cancer$dead == 1, 1, 0)
cancer <- cancer[cancer$time >= 365 | cancer$dead == 1, ]

# Estimate all 7 contrasts in one call
result <- aipw_estimate(
  data       = cancer,
  outcome    = "dead_365",
  treatment  = "ecog_bin",
  covariates = c("age", "sex"),
  sl_library_PS = c("SL.glm", "SL.earth"),
  sl_library_OR = c("SL.glm", "SL.earth")
)
result
```

```
=== AIPW Causal Estimates ===

 measure estimate   lower   upper
      RD   0.1509 -0.0180  0.3198
     ERR   0.2794 -0.0515  0.7259
      RR   1.2794  0.9485  1.7259
      OR   1.9042  0.9404  3.8557
     NNT   6.6267  2.1636 20.2957
      SR   0.6719  0.4447  1.0152
      VE  -0.2794 -0.7259  0.0515
   E[Y1]   0.6919  0.5869  0.7969
   E[Y0]   0.5410  0.4463  0.6357
```

## Key Functions

| Function | Description |
|---|---|
| `aipw_estimate()` | Cross-fitted AIPW estimator with influence-function-based CIs |
| `aipw_bootstrap()` | Bootstrap confidence intervals as an alternative to EIF-based inference |
| `plot_estimates()` | Forest-style plots for each contrast |
| `plot_covariate_diagnostics()` | Density plots to assess overlap and covariate balance |
| `available_contrasts()` | List all supported causal contrasts |
| `available_sl_libraries()` | List all supported SuperLearner algorithms |

## Features

### Separate SuperLearner libraries

Specify different algorithms for the propensity score and outcome models:

```r
aipw_estimate(
  ...,
  sl_library_PS = c("SL.glm", "SL.ranger"),     # for P(A=1|X)
  sl_library_OR = c("SL.xgboost", "SL.earth")   # for E[Y|A,X]
)
```

### Select specific contrasts

```r
# Only report risk difference and risk ratio
aipw_estimate(..., contrast = c("RD", "RR"))
```

### Diagnostics

```r
# Check overlap between treatment groups
plot_covariate_diagnostics(
  data = cancer,
  treatment = "ecog_bin",
  covariates = c("age", "sex"),
  labels = c("Asymptomatic", "Has impairment")
)
```

### Visualization

```r
result <- aipw_estimate(...)

# Plot all measures
plot_estimates(result)

# Plot specific measures and save
plot_estimates(result, measures = c("RD", "NNT"), save_dir = "figures/")
```

## Getting Started

The best place to get started is the Cancer Tutorial vignette, which walks through a complete analysis of the lung cancer dataset from the `survival` package:

```r
vignette("Cancer_tutorial", package = "MultiContrasts")
```

## Issues

If you encounter any bugs or have any specific feature requests, please [file an issue](https://github.com/zhiyuemo/MultiContrasts/issues).

## Contributions

Contributions are very welcome. Interested contributors should fork the repository and submit a pull request. Please ensure that `devtools::check()` passes with 0 errors and 0 warnings before submitting.

## Citation

After using the `MultiContrasts` R package, please cite the following:

```bibtex
@software{mo2026multicontrasts,
  author = {Mo, Zhiyue},
  title = {{MultiContrasts}: Doubly Robust Estimation of Multiple Causal Contrasts},
  year = {2026},
  howpublished = {\url{https://github.com/zhiyuemo/MultiContrasts}},
  note = {{R} package version 0.0.0.9000}
}
```

## License

© 2026 Zhiyue Mo

The contents of this repository are distributed under the MIT license. See file [LICENSE](LICENSE) for details.

## References

- Robins JM, Rotnitzky A, Zhao LP. Estimation of regression coefficients when some regressors are not always observed. *Journal of the American Statistical Association*. 1994;89(427):846-866.
- van der Laan MJ, Polley EC, Hubbard AE. Super Learner. *Statistical Applications in Genetics and Molecular Biology*. 2007;6(1).
- Chernozhukov V, Chetverikov D, Demirer M, et al. Double/debiased machine learning for treatment and structural parameters. *The Econometrics Journal*. 2018;21(1):C1-C68.
- Boughdiri N, et al. A unified framework for transporting causal measures under covariate shift. 2025.
