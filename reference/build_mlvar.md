# Build a Multilevel Vector Autoregression (mlVAR) network

Estimates three networks from ESM/EMA panel data, matching
[`mlVAR::mlVAR()`](https://rdrr.io/pkg/mlVAR/man/mlVAR.html) with
`estimator = "lmer"`, `temporal = "fixed"`, `contemporaneous = "fixed"`
at machine precision: (1) a directed temporal network of fixed-effect
lagged regression coefficients, (2) an undirected contemporaneous
network of partial correlations among residuals, and (3) an undirected
between-subjects network of partial correlations derived from the
person-mean fixed effects.

\#' @details The algorithm follows mlVAR's lmer pipeline exactly:

1.  Drop rows with NA in id/day/beep and optionally grand-mean
    standardize each variable.

2.  Expand the per-(id, day) beep grid and right-join original values,
    producing the augmented panel (`augData`).

3.  Add within-person lagged predictors (`L1_*`) and person-mean
    predictors (`PM_*`).

4.  For each outcome variable fit
    `lmer(y ~ within + between-except-own-PM + (1 | id))` with
    `REML = FALSE`. Collect the fixed-effect temporal matrix `B`,
    between-effect matrix `Gamma`, random-intercept SDs (`mu_SD`), and
    lmer residual SDs.

5.  Contemporaneous network:
    `cor2pcor(D %*% cov2cor(cor(resid)) %*% D)`.

6.  Between-subjects network:
    `cor2pcor(pseudoinverse(forcePositive(D (I - Gamma))))`.

Validated to machine precision (max_diff \< 1e-10) against
[`mlVAR::mlVAR()`](https://rdrr.io/pkg/mlVAR/man/mlVAR.html) on 25 real
ESM datasets from `openesm` and 20 simulated configurations (seeds
201-220). See `tmp/mlvar_equivalence_real20.R` and
`tmp/mlvar_equivalence_20seeds.R`.

## Usage

``` r
build_mlvar(
  data,
  vars,
  id,
  day = NULL,
  beep = NULL,
  lag = 1L,
  standardize = FALSE
)
```

## Arguments

- data:

  A `data.frame` containing the panel data.

- vars:

  Character vector of variable column names to model.

- id:

  Character string naming the person-ID column.

- day:

  Character string naming the day/session column, or `NULL`. When
  provided, lag pairs are only formed within the same day.

- beep:

  Character string naming the measurement-occasion column, or `NULL`.
  When `NULL`, row position within each (id, day) is used.

- lag:

  Integer. The lag order (default 1).

- standardize:

  Logical. If `TRUE`, each variable is grand-mean centered and divided
  by its pooled SD *before* augmentation. Default `FALSE`, matching
  `mlVAR::mlVAR(scale = FALSE)` — the only setting for which numerical
  equivalence has been validated.

## Value

A dual-class `c("net_mlvar", "netobject_group")` object — a named list
of three full netobjects, one per network, plus model-level metadata
stored as attributes. Each element is a standard
`c("netobject", "cograph_network")` weight-matrix wrapper (no raw
`$data`), so [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coefs()`](https://saqr.me/Nestimate/reference/coefs.md), and
`cograph::splot(fit$temporal)` work directly. The three constituents are
matrix-wrapped and carry no underlying panel data, so data-resampling
verbs such as
[`bootstrap_network()`](https://saqr.me/Nestimate/reference/bootstrap_network.md)
(and reliability/stability) cannot iterate over them — extract a single
constituent and rebuild via
[`build_network()`](https://saqr.me/Nestimate/reference/build_network.md)
if you need those. Structure:

- `fit$temporal`:

  Directed netobject for the `d x d` matrix of fixed-effect lagged
  coefficients. `$weights[i, j]` is the effect of variable j at t-lag on
  variable i at t. `method = "mlvar_temporal"`, `directed = TRUE`.

- `fit$contemporaneous`:

  Undirected netobject for the `d x d` partial-correlation network of
  within-person lmer residuals. `method = "mlvar_contemporaneous"`,
  `directed = FALSE`.

- `fit$between`:

  Undirected netobject for the `d x d` partial-correlation network of
  person means, derived from `D (I - Gamma)`.
  `method = "mlvar_between"`, `directed = FALSE`.

- `attr(fit, "coefs")` /
  [`coefs()`](https://saqr.me/Nestimate/reference/coefs.md):

  Tidy `data.frame` with one row per `(outcome, predictor)` pair and
  columns `outcome`, `predictor`, `beta`, `se`, `t`, `p`, `ci_lower`,
  `ci_upper`, `significant`. Filter, sort, or plot with base R or the
  tidyverse. Retrieve with `coefs(fit)`.

- `attr(fit, "n_obs")`:

  Number of rows in the augmented panel after na.omit.

- `attr(fit, "n_subjects")`:

  Number of unique subjects remaining.

- `attr(fit, "lag")`:

  Lag order used.

- `attr(fit, "standardize")`:

  Logical; whether pre-augmentation standardization was applied.

## See also

[`build_network()`](https://saqr.me/Nestimate/reference/build_network.md)

## Examples

``` r
if (FALSE) { # \dontrun{
d <- simulate_data("mlvar", seed = 1)
fit <- build_mlvar(d, vars = attr(d, "vars"),
                   id = "id", day = "day", beep = "beep")
print(fit)
summary(fit)
} # }
```
