# Cluster Sequences by Dissimilarity

Clusters wide-format sequences using pairwise string dissimilarity and
either PAM (Partitioning Around Medoids) or hierarchical clustering.
Supports 9 distance metrics including temporal weighting for Hamming
distance. When the stringdist package is available, uses C-level
distance computation for 100-1000x speedup on edit distances.

## Usage

``` r
build_clusters(
  data,
  k,
  dissimilarity = "hamming",
  method = "pam",
  na_syms = c("*", "%"),
  weighted = FALSE,
  lambda = 1,
  seed = NULL,
  q = 2L,
  p = 0.1,
  covariates = NULL,
  estimator = c("firth", "multinom", "chisq"),
  ...
)
```

## Arguments

- data:

  Input data. Accepts multiple formats:

  data.frame / matrix

  :   Wide-format sequences (rows = sequences, columns = time points,
      values = state names).

  netobject

  :   A network object from
      [`build_network`](https://saqr.me/Nestimate/reference/build_network.md).
      Extracts the stored sequence data. Only valid for sequence-based
      methods (relative, frequency, co_occurrence, attention).

  tna

  :   A tna model from the tna package. Decodes the integer-encoded
      sequence data using stored labels.

  cograph_network

  :   A cograph network object. Extracts the stored sequence data.

- k:

  Integer. Number of clusters (must be between 2 and `nrow(data) - 1`).

- dissimilarity:

  Character. Distance metric. One of `"hamming"`, `"osa"` (optimal
  string alignment), `"lv"` (Levenshtein), `"dl"` (Damerau-Levenshtein),
  `"lcs"` (longest common subsequence), `"qgram"`, `"cosine"`,
  `"jaccard"`, `"jw"` (Jaro-Winkler). Default: `"hamming"`.

- method:

  Character. Clustering method. `"pam"` for Partitioning Around Medoids,
  or a hierarchical method: `"ward.D2"`, `"ward.D"`, `"complete"`,
  `"average"`, `"single"`, `"mcquitty"`, `"median"`, `"centroid"`.
  Default: `"pam"`.

- na_syms:

  Character vector. Symbols treated as missing values. Default:
  `c("*", "%")`.

  **Missing-value distance rule:** after symbols are converted to `NA`,
  missing values are encoded as a single comparable sentinel state –
  *not* pairwise-deleted. Two missing values in the same position match
  (distance contribution 0); a missing value paired with any observed
  state mismatches (distance contribution 1 for Hamming, etc.). This is
  the conventional behaviour for aligned sequence matrices because
  pairwise deletion would change the effective length of every pair and
  break the metric. If you want pairwise deletion or a different
  missing-value semantic, drop or recode the missing cells before
  passing the data in.

- weighted:

  Logical. Apply exponential decay weighting to Hamming distance
  positions? Only valid when `dissimilarity = "hamming"`. Default:
  `FALSE`.

- lambda:

  Numeric. Non-negative decay rate for weighted Hamming. Higher values
  weight earlier positions more strongly. Default: 1.

- seed:

  Integer or NULL. Random seed for reproducibility. Default: `NULL`.

- q:

  Integer. Size of q-grams for `"qgram"`, `"cosine"`, and `"jaccard"`
  distances. Default: `2L`.

- p:

  Numeric. Winkler prefix penalty for Jaro-Winkler distance. Must be
  between 0 and 0.25. Default: `0.1`.

- covariates:

  Optional. Post-hoc covariate analysis of cluster membership. Accepts:

  string

  :   Single column name, e.g. `"Age"`. Resolved against `x$metadata`
      (and `x$data`) for `netobject` or `cograph_network` input.

  character vector

  :   `c("Age", "Gender")`, same lookup.

  formula

  :   `~ Age + Gender`, same lookup; supports `"Age + Gender"` string
      form too.

  data.frame

  :   All columns used as covariates verbatim; must have one row per
      sequence.

  NULL

  :   No covariate analysis (default).

  For `netobject` or `cograph_network` input, names are resolved against
  `$metadata` first and then non-state columns of `$data`, so a typical
  call looks like
  `build_clusters(net, k = 3, covariates = "session_label")` without
  pre-extracting a data.frame. `tna` input requires the data.frame form.
  Results are stored in `$covariates`.

- estimator:

  Multinomial logit fitter for the covariate analysis. `"firth"`
  (default) uses Firth's penalised likelihood via
  [`brglm2::brmultinom`](https://rdrr.io/pkg/brglm2/man/brmultinom.html)
  — bias-reduced and finite under quasi-complete separation.
  `"multinom"` uses classical ML via
  [`nnet::multinom`](https://rdrr.io/pkg/nnet/man/multinom.html); emits
  a warning because rare-cell separation produces astronomical ORs with
  degenerate CIs (silent failure). `"chisq"` runs WeightedCluster-style
  descriptive tests (chi-square + Cramer's V + standardized adjusted
  residuals for factors; Kruskal-Wallis + eta-squared for numerics).

- ...:

  Unsupported. Supplying unused arguments raises an error.

## Value

An object of class `"net_clustering"` containing:

- data:

  The original input data.

- k:

  Number of clusters.

- assignments:

  Named integer vector of cluster assignments.

- silhouette:

  Overall average silhouette width.

- sizes:

  Named integer vector of cluster sizes.

- method:

  Clustering method used.

- dissimilarity:

  Distance metric used.

- distance:

  The computed dissimilarity matrix (`dist` object).

- medoids:

  Integer vector of medoid row indices (PAM only; NULL for hierarchical
  methods).

- seed:

  Seed used (or NULL).

- weighted:

  Logical, whether weighted Hamming was used.

- lambda:

  Lambda value used (0 if not weighted).

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","C","A","B"), V2 = c("B","C","A","B","A"),
                   V3 = c("C","A","B","C","B"))
cl <- build_clusters(seqs, k = 2)
cl
#> Sequence Clustering [pam]
#>   Sequences: 5  |  Clusters: 2
#>   Dissimilarity: hamming
#>   Quality: silhouette = 0.600
#> 
#>   Cluster  N          Mean within-dist  Medoid
#>   1        2 (40.0%)  0.000             4
#>   2        3 (60.0%)  2.000             5
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:3], 20, TRUE), V2 = sample(LETTERS[1:3], 20, TRUE),
  V3 = sample(LETTERS[1:3], 20, TRUE), V4 = sample(LETTERS[1:3], 20, TRUE)
)
cl <- build_clusters(seqs, k = 2)
print(cl)
#> Sequence Clustering [pam]
#>   Sequences: 20  |  Clusters: 2
#>   Dissimilarity: hamming
#>   Quality: silhouette = 0.266
#> 
#>   Cluster  N           Mean within-dist  Medoid
#>   1        9 (45.0%)   2.389             13
#>   2        11 (55.0%)  2.036             20
summary(cl)
#> Sequence Clustering Summary
#>   Method:        pam 
#>   Dissimilarity: hamming 
#>   Silhouette:    0.2659 
#> 
#> Per-cluster statistics:
#>  cluster size mean_within_dist
#>        1    9         2.388889
#>        2   11         2.036364
#>   cluster size mean_within_dist
#> 1       1    9         2.388889
#> 2       2   11         2.036364
# }
```
