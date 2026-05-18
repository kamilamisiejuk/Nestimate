# Compute Centrality Measures for a Network

Computes centrality measures from a `netobject`, `netobject_group`, or
`cograph_network`. For directed networks the default measures are
InStrength, OutStrength, and Betweenness. For undirected networks the
defaults are Closeness and Betweenness.

## Usage

``` r
net_centrality(x, measures = NULL, loops = FALSE, centrality_fn = NULL, ...)
```

## Arguments

- x:

  A `netobject`, `netobject_group`, or `cograph_network`.

- measures:

  Character vector. Centrality measures to compute. Built-in:
  `"InStrength"`, `"OutStrength"`, `"Betweenness"`, `"InCloseness"`,
  `"OutCloseness"`, `"Closeness"`. `"Closeness"` is defined only for
  undirected networks; `"InCloseness"`/`"OutCloseness"` only for
  directed networks (requesting the wrong one for the network's
  directedness is an error). Default depends on directedness.

- loops:

  Logical. Include self-loops (diagonal) in computation? Default:
  `FALSE`.

- centrality_fn:

  Optional function. Custom centrality function that takes a weight
  matrix and returns a named list of centrality vectors.

- ...:

  Additional arguments (ignored).

## Value

For a `netobject`: a data frame with node names as rows and centrality
measures as columns. For a `netobject_group`: a named list of such data
frames (one per group).

## Examples

``` r
seqs <- data.frame(
  V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
  V3 = c("C","A","C","B"))
net <- build_network(seqs, method = "relative")
net_centrality(net)
#> centralities computed excluding loops (diagonal). Pass `loops = TRUE` to include self-transitions.
#>   InStrength OutStrength Betweenness
#> A          1           1         0.5
#> B          1           1         0.5
#> C          1           1         0.5
```
