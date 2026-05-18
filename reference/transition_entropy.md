# Transition Entropy of a Markov Chain

Computes per-state branching entropy, stationary entropy, and the
chain-level entropy rate of a Markov transition process. The entropy
rate is the Shannon-McMillan-Breiman per-step uncertainty of
trajectories under the stationary distribution; it is the canonical
information-theoretic summary of a transition matrix.

## Usage

``` r
transition_entropy(x, base = 2, normalize = TRUE)
```

## Arguments

- x:

  A `netobject`, `cograph_network`, `tna` object, row-stochastic numeric
  transition matrix, or a wide sequence data.frame (rows = actors,
  columns = time-steps; a relative transition network is built
  automatically). Group dispatch on `netobject_group`.

- base:

  Numeric. Logarithm base. `2` (default) for bits, `exp(1)` for nats,
  `10` for hartleys.

- normalize:

  Logical. If `TRUE` (default), rows that do not sum to 1 are normalised
  automatically (with a warning).

## Value

An object of class `"net_transition_entropy"` with:

- row_entropy:

  Named numeric vector, length \\n\\. Per-state branching entropy
  \\H(P\_{i\cdot}) = -\sum_j P\_{ij} \log P\_{ij}\\.

- row_entropy_norm:

  Named numeric vector. `row_entropy` divided by the ceiling \\\log_b
  n\\ (in \\\[0, 1\]\\; all zeros when \\n = 1\\).

- stationary:

  Named numeric vector. Stationary distribution \\\pi\\.

- stationary_entropy:

  Scalar. \\H(\pi) = -\sum_i \pi_i \log \pi_i\\ - the entropy of \\\pi\\
  treated as an i.i.d. distribution. Upper bound on the entropy rate.

- stationary_entropy_norm:

  Scalar. `stationary_entropy` divided by the ceiling \\\log_b n\\.

- entropy_rate:

  Scalar. \\h(P) = \sum_i \pi_i H(P\_{i\cdot})\\ - the
  Shannon-McMillan-Breiman entropy rate.

- entropy_rate_norm:

  Scalar. `entropy_rate` divided by the ceiling \\\log_b n\\.

- redundancy:

  Scalar. \\H(\pi) - h(P)\\, the entropy deficit attributable to serial
  dependence; zero for an i.i.d. chain (rows of \\P\\ all equal
  \\\pi\\).

- redundancy_norm:

  Scalar. The *relative* redundancy \\(H(\pi) - h(P)) / H(\pi)\\ (the
  fraction of the stationary entropy removed by order-1 memory), **not**
  `redundancy` divided by \\\log_b n\\; `0` when \\H(\pi) = 0\\.

- max_entropy:

  Scalar. The normalising ceiling \\\log_b n\\.

- base:

  Logarithm base used.

- states:

  Character vector of state names.

## Details

Convention \\0 \log 0 := 0\\ is applied, so absorbing or deterministic
rows contribute zero per-row entropy. The chain need not be irreducible;
\\\pi\\ is computed from the eigendecomposition of \\P^\top\\ as
elsewhere in the package. For non-ergodic chains the returned \\\pi\\ is
one stationary distribution among many - interpret with the help of
[`chain_structure`](https://saqr.me/Nestimate/reference/chain_structure.md).

The relation \\h(P) \leq H(\pi)\\ holds with equality iff successive
states are independent. The deficit \\H(\pi) - h(P)\\ is reported as
`redundancy` - a measure of how much memory the chain has at order 1.

## References

Cover, T.M. & Thomas, J.A. (2006). *Elements of Information Theory*, 2nd
ed., chapter 4. Wiley.

Shannon, C.E. (1948). A mathematical theory of communication. *Bell
System Technical Journal*, 27, 379-423.

## See also

[`markov_stability`](https://saqr.me/Nestimate/reference/markov_stability.md),
[`passage_time`](https://saqr.me/Nestimate/reference/passage_time.md),
[`markov_order_test`](https://saqr.me/Nestimate/reference/markov_order_test.md),
[`chain_structure`](https://saqr.me/Nestimate/reference/chain_structure.md)

## Examples

``` r
# \donttest{
net <- build_network(as.data.frame(trajectories), method = "relative")
te  <- transition_entropy(net)
print(te)
#> Transition Entropy (3 states, bits; ceiling = 1.585)
#> 
#>                           raw            normalised
#>   Entropy rate    h(P)  = 1.247 bits    0.787
#>   Stationary    H(pi)  = 1.501 bits    0.947
#>   Redundancy   H(pi)-h = 0.255 bits    0.170
#> 
#> Normalised: h(P) and H(pi) are raw / log_2(n_states) (0 = deterministic, 1 = uniform);
#>   redundancy is the relative redundancy (H(pi) - h(P)) / H(pi), not raw / log_2(n_states).
summary(te)
#> Transition Entropy Summary (bits)
#> 
#> Per-state contribution to h(P):
#>       state stationary row_entropy row_entropy_norm contribution
#>     Average      0.443       1.354            0.855        0.600
#>      Active      0.372       1.040            0.656        0.387
#>  Disengaged      0.185       1.403            0.885        0.260
#>  contribution_pct
#>              48.1
#>              31.0
#>              20.8
#> 
#> Chain-level summary:
#>               quantity   raw normalised
#>      entropy_rate h(P) 1.247      0.787
#>       stationary H(pi) 1.501      0.947
#>  redundancy H(pi)-h(P) 0.255      0.170
#>       ceiling log_2(n) 1.585      1.000
#> 
#> Normalised: row entropy, h(P) and H(pi) are raw / log_2(n_states), in [0, 1];
#>   redundancy is the relative redundancy (H(pi) - h(P)) / H(pi).
plot(te)

# }
```
