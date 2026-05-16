# Persistent Homology

Builds clique complexes at decreasing inclusive weight thresholds and
records the Betti curve. The returned persistence table is a heuristic
pairing derived from changes in Betti counts across the threshold grid;
it is not a full boundary-matrix persistent homology decomposition of
individual homology classes.

## Usage

``` r
persistent_homology(x, n_steps = 20L, max_dim = 3L)
```

## Arguments

- x:

  A square matrix, `tna`, or `netobject`.

- n_steps:

  Number of filtration steps (default 20).

- max_dim:

  Maximum simplex dimension to track (default 3).

## Value

A `persistent_homology` object with:

- betti_curve:

  Data frame: `threshold`, `dimension`, `betti` at each filtration step.

- persistence:

  Data frame of birth-death pairs: `dimension`, `birth`, `death`,
  `persistence`. These intervals summarize Betti-count changes across
  the sampled thresholds.

- thresholds:

  Numeric vector of filtration thresholds.

## Examples

``` r
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
ph <- persistent_homology(mat, n_steps = 10)
print(ph)
#> Persistent Homology
#>   10 filtration steps [0.6000 → 0.0060]
#>   Features: b0: 2 (1 persistent) 
#>   Longest-lived:
#>     b0: 0.6000 → 0.0000 (life: 0.6000)
#>     b0: 0.6000 → 0.4680 (life: 0.1320)
```
