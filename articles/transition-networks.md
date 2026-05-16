# Network Estimation and Analysis with Nestimate

`Nestimate` is a unified framework for estimating, validating, and
comparing networks from sequential and cross-sectional data. It
implements two complementary paradigms: **Transition Network Analysis
(TNA)**, which models the relational dynamics of temporal processes as
weighted directed networks using stochastic Markov models; and
**Psychological Network Analysis (PNA)**, which estimates the
conditional dependency structure among variables using regularized
partial correlations and graphical models. Both paradigms share the same
[`build_network()`](https://saqr.me/Nestimate/reference/build_network.md)
interface, the same validation engine (bootstrap, permutation,
centrality stability), and the same output format — enabling researchers
to apply a consistent analytic workflow across fundamentally different
data types.

This vignette demonstrates both paradigms, covering network estimation,
statistical validation, data-driven clustering, and group comparison.

## Part I: Transition Network Analysis

### Theoretical Grounding

TNA uses stochastic process modeling to capture the dynamics of temporal
processes via Markov models. Markov models align with the view that a
temporal process is an outcome of a stochastic data generating process
that produces various network configurations or patterns based on rules,
constraints, or guiding principles. The transitions are governed by a
stochastic process: the specific ways in which the system changes or
evolves is rather random and therefore cannot be strictly determined.
That is, the transitions are probabilistically dependent on preceding
states — a method that assumes events are probabilistically dependent on
the preceding ones.

The main principle of TNA is representing the transition matrix between
events as a graph to take full advantage of graph theory potentials and
the wealth of network analysis. TNA brings network measures at the node,
edge, and graph level; pattern mining through dyads, triads, and
communities; clustering of sub-networks into typical behavioral
strategies; and rigorous statistical validation at each edge through
bootstrapping, permutation, and case-dropping techniques. Such
statistical rigor that brings validation and hypothesis testing at each
step of the analysis offers a method for researchers to build, verify,
and advance existing theories on the basis of a robust scientific
approach.

For a comprehensive introduction to TNA in learning analytics, see Saqr
et al. (2025a) and Saqr (2024). For detailed tutorials on TNA clustering
and heterogeneity analysis, see López-Pernas et al. (2025).

### Data

The `human_long` dataset contains 10,796 coded human interaction turns
from 429 human-AI pair programming sessions across 34 projects. Each row
is a single turn, with `code` recording the interaction type,
`session_id` identifying the session, and `timestamp` providing temporal
ordering. For a detailed description of the dataset and coding scheme,
see [Saqr
(2026)](https://saqr.me/blog/2026/human-ai-interaction-cograph/).

``` r

library(Nestimate)

# Subsample for vignette speed (CRAN build-time limit)
set.seed(1)
keep <- sample(unique(human_long$session_id), 100)
human_sub <- human_long[human_long$session_id %in% keep, ]

head(human_sub)
#>     message_id   project   session_id  timestamp session_date      code
#> 395       2902 Project_7 0605767ae57f 1772229600   2026-02-28   Specify
#> 396       2902 Project_7 0605767ae57f 1772229600   2026-02-28   Command
#> 397       2902 Project_7 0605767ae57f 1772229600   2026-02-28   Request
#> 398       2902 Project_7 0605767ae57f 1772229600   2026-02-28   Specify
#> 399       2903 Project_7 0605767ae57f 1772229600   2026-02-28 Interrupt
#> 400       2905 Project_7 0605767ae57f 1772229600   2026-02-28   Command
#>           cluster code_order order_in_session
#> 395     Directive          1                1
#> 396     Directive          2                2
#> 397     Directive          3                3
#> 398     Directive          4                4
#> 399 Metacognitive          1                5
#> 400     Directive          1                8
```

The dataset is in long format: `code` records what happened,
`session_id` who did it, and `timestamp` when. Additional columns like
`project` and `cluster` are automatically preserved as metadata for
downstream covariate analysis.

### Building Networks

Building networks in Nestimate is a single step:
[`build_network()`](https://saqr.me/Nestimate/reference/build_network.md)
is the universal entry point for all network estimation. It accepts
long-format event data directly with three key parameters:

- **`action`**: the column containing state labels
- **`actor`**: the column identifying sequences (one sequence per actor)
- **`time`**: the column providing temporal ordering

Under the hood,
[`build_network()`](https://saqr.me/Nestimate/reference/build_network.md)
automatically converts the long-format event log into wide-format
sequences, handling chronological ordering, session detection, and
metadata preservation. You can also call
[`prepare()`](https://saqr.me/Nestimate/reference/prepare.md) directly
to inspect or reuse the processed data before passing it to
[`build_network()`](https://saqr.me/Nestimate/reference/build_network.md).

#### Transition Network (TNA)

The standard TNA method estimates a first-order Markov model from
sequence data. Given a sequence of events, the transition probability
$`P(v_j | v_i)`$ is estimated as the ratio of observed transitions from
state $`v_i`$ to state $`v_j`$ to the total number of outgoing
transitions from $`v_i`$. These probabilities are assembled into a
**transition matrix** $`T`$, where each element $`T_{ij}`$ represents
the estimated probability of transitioning from $`v_i`$ to $`v_j`$ (Saqr
et al., 2025a).

``` r

net_tna <- build_network(human_sub, method = "tna",
                         action = "code", actor = "session_id",
                         time = "timestamp")
print(net_tna)
#> Transition Network (relative probabilities) [directed]
#>   Weights: [0.019, 0.577]  |  mean: 0.111
#> 
#>   Weight matrix:
#>             Command Correct Frustrate Inquire Interrupt Refine Request Specify
#>   Command     0.209   0.091     0.040   0.067     0.037  0.047   0.158   0.295
#>   Correct     0.099   0.109     0.134   0.054     0.059  0.109   0.124   0.282
#>   Frustrate   0.078   0.150     0.204   0.083     0.019  0.209   0.107   0.117
#>   Inquire     0.155   0.155     0.140   0.166     0.109  0.041   0.093   0.098
#>   Interrupt   0.298   0.055     0.055   0.193     0.144  0.039   0.066   0.077
#>   Refine      0.042   0.090     0.090   0.066     0.048  0.096   0.114   0.416
#>   Request     0.102   0.019     0.070   0.074     0.051  0.033   0.033   0.577
#>   Specify     0.291   0.069     0.076   0.067     0.176  0.052   0.064   0.162
#>   Verify      0.230   0.100     0.110   0.100     0.050  0.140   0.080   0.110
#>             Verify
#>   Command    0.056
#>   Correct    0.030
#>   Frustrate  0.034
#>   Inquire    0.041
#>   Interrupt  0.072
#>   Refine     0.036
#>   Request    0.042
#>   Specify    0.043
#>   Verify     0.080 
#> 
#>   Initial probabilities:
#>   Specify       0.664  ████████████████████████████████████████
#>   Command       0.216  █████████████
#>   Request       0.048  ███
#>   Correct       0.040  ██
#>   Interrupt     0.016  █
#>   Frustrate     0.008  
#>   Refine        0.008  
#>   Inquire       0.000  
#>   Verify        0.000
```

#### Frequency Network (FTNA)

The frequency method preserves raw transition counts rather than
normalizing to conditional probabilities. This is useful when absolute
frequencies matter — a transition occurring 500 times from a common
state may be more practically important than one occurring 5 times from
a rare state, even if the latter has a higher conditional probability.

``` r

net_ftna <- build_network(human_sub, method = "ftna",
                          action = "code", actor = "session_id",
                          time = "timestamp")
print(net_ftna)
#> Transition Network (frequency counts) [directed]
#>   Weights: [4.000, 169.000]  |  mean: 28.062
#> 
#>   Weight matrix:
#>             Command Correct Frustrate Inquire Interrupt Refine Request Specify
#>   Command        90      39        17      29        16     20      68     127
#>   Correct        20      22        27      11        12     22      25      57
#>   Frustrate      16      31        42      17         4     43      22      24
#>   Inquire        30      30        27      32        21      8      18      19
#>   Interrupt      54      10        10      35        26      7      12      14
#>   Refine          7      15        15      11         8     16      19      69
#>   Request        22       4        15      16        11      7       7     124
#>   Specify       169      40        44      39       102     30      37      94
#>   Verify         23      10        11      10         5     14       8      11
#>             Verify
#>   Command       24
#>   Correct        6
#>   Frustrate      7
#>   Inquire        8
#>   Interrupt     13
#>   Refine         6
#>   Request        9
#>   Specify       25
#>   Verify         8 
#> 
#>   Initial probabilities:
#>   Specify       0.664  ████████████████████████████████████████
#>   Command       0.216  █████████████
#>   Request       0.048  ███
#>   Correct       0.040  ██
#>   Interrupt     0.016  █
#>   Frustrate     0.008  
#>   Refine        0.008  
#>   Inquire       0.000  
#>   Verify        0.000
```

#### Attention Network (ATNA)

The attention method applies temporal decay weighting, giving more
importance to recent transitions within each sequence. The `lambda`
parameter controls the decay rate: higher values produce faster decay.
This captures the idea that later events in a process may be more
indicative of the underlying dynamics than early ones.

``` r

net_atna <- build_network(human_sub, method = "atna",
                          action = "code", actor = "session_id",
                          time = "timestamp")
print(net_atna)
#> Attention Network (decay-weighted transitions) [directed]
#>   Weights: [2.719, 81.786]  |  mean: 15.834
#> 
#>   Weight matrix:
#>             Command Correct Frustrate Inquire Interrupt Refine Request Specify
#>   Command    50.398  21.700    12.450  16.556    15.164 11.968  31.860  68.371
#>   Correct    13.595  13.498    14.098   8.029     7.025 12.674  11.888  30.257
#>   Frustrate  10.506  15.087    23.142   9.736     3.927 20.308  12.254  18.614
#>   Inquire    16.073  17.002    13.903  16.691    11.297  4.845   9.755  14.765
#>   Interrupt  25.913   6.000     5.844  17.104    17.124  3.850   6.353  12.338
#>   Refine      6.765   9.389     9.894   6.244     4.124 10.875  10.187  33.191
#>   Request    14.348   4.176     8.718  10.071    10.654  6.075   6.291  53.456
#>   Specify    81.786  22.700    26.216  25.113    44.442 18.861  28.052  65.762
#>   Verify     11.587   5.170     5.991   5.540     2.719  7.318   5.262   8.894
#>             Verify
#>   Command   13.721
#>   Correct    3.554
#>   Frustrate  3.574
#>   Inquire    4.062
#>   Interrupt  6.741
#>   Refine     3.204
#>   Request    5.989
#>   Specify   15.915
#>   Verify     3.943 
#> 
#>   Initial probabilities:
#>   Specify       0.664  ████████████████████████████████████████
#>   Command       0.216  █████████████
#>   Request       0.048  ███
#>   Correct       0.040  ██
#>   Interrupt     0.016  █
#>   Frustrate     0.008  
#>   Refine        0.008  
#>   Inquire       0.000  
#>   Verify        0.000
```

#### Co-occurrence Network from Binary Data

When the data is binary (0/1) — as is common in learning analytics where
activities are coded as present or absent within time windows —
[`build_network()`](https://saqr.me/Nestimate/reference/build_network.md)
automatically detects the format and uses co-occurrence analysis. The
resulting undirected network captures which events tend to co-occur.

``` r

data(learning_activities)
net <- build_network(learning_activities, method = "cna", actor = "student")
print(net)
#> Co-occurrence Network [undirected]
#>   Weights: [2681.000, 6262.000]  |  mean: 3867.381
#> 
#>   Weight matrix:
#>           Reading Video Forum Quiz Coding Review
#>   Reading    6100  3169  3211 2891   3138   3113
#>   Video      3169  6262  3183 2903   3278   3290
#>   Forum      3211  3183  5692 2942   2958   3045
#>   Quiz       2891  2903  2942 5379   2681   2725
#>   Coding     3138  3278  2958 2681   5886   3183
#>   Review     3113  3290  3045 2725   3183   6186
```

#### Window-based TNA (WTNA)

The [`wtna()`](https://saqr.me/Nestimate/reference/wtna.md) function
computes networks from one-hot encoded (binary) data using temporal
windowing. It supports three modes:

- **`"transition"`**: directed transitions between consecutive windows
- **`"cooccurrence"`**: undirected co-occurrence within windows
- **`"both"`**: a mixed network combining transitions and co-occurrences

``` r

net_wtna <- wtna(learning_activities, actor = "student",
                 method = "transition", type = "frequency")
print(net_wtna)
#> Network (method: wtna_transition) [directed]
#>   Weights: [2493.000, 4314.000]  |  mean: 3025.056
#> 
#>   Weight matrix:
#>           Reading Video Forum Quiz Coding Review
#>   Reading    4297  2903  2908 2638   2935   3008
#>   Video      2921  4314  2943 2651   2909   3018
#>   Forum      2863  2875  3881 2675   2747   2902
#>   Quiz       2691  2654  2706 3695   2518   2542
#>   Coding     2814  3094  2767 2559   4005   2953
#>   Review     2971  3113  2824 2493   2909   4206 
#> 
#>   Initial probabilities:
#>   Reading       0.209  ████████████████████████████████████████
#>   Coding        0.194  █████████████████████████████████████
#>   Video         0.179  ██████████████████████████████████
#>   Review        0.151  █████████████████████████████
#>   Quiz          0.141  ███████████████████████████
#>   Forum         0.126  ████████████████████████
```

``` r

net_wtna_rel <- wtna(learning_activities, method = "transition", type = "relative")
print(net_wtna_rel)
#> Network (method: wtna_transition) [directed]
#>   Weights: [0.136, 0.225]  |  mean: 0.167
#> 
#>   Weight matrix:
#>           Reading Video Forum  Quiz Coding Review
#>   Reading   0.225 0.157 0.155 0.142  0.159  0.162
#>   Video     0.156 0.225 0.157 0.142  0.159  0.161
#>   Forum     0.160 0.161 0.210 0.150  0.156  0.162
#>   Quiz      0.161 0.160 0.161 0.213  0.151  0.153
#>   Coding    0.157 0.169 0.152 0.141  0.220  0.162
#>   Review    0.162 0.167 0.152 0.136  0.159  0.224 
#> 
#>   Initial probabilities:
#>   Reading       0.500  ████████████████████████████████████████
#>   Quiz          0.500  ████████████████████████████████████████
#>   Video         0.000  
#>   Forum         0.000  
#>   Coding        0.000  
#>   Review        0.000
```

#### Mixed Network (Transitions + Co-occurrences)

Since states can co-occur within the same window *and* follow each other
across windows, a mixed network captures both relationships
simultaneously — modeling co-activity and temporal succession in a
single structure.

``` r

net_wtna_mixed <- wtna(learning_activities, method = "both", type = "relative")
print(net_wtna_mixed)
#> Mixed Window TNA (transition + co-occurrence)
#> -- Transition (directed) --
#>   Nodes: 6  |  Edges: 36
#> -- Co-occurrence (undirected) --
#>   Nodes: 6  |  Edges: 21
```

### Validation

Most research on networks or process mining uses descriptive methods.
The validation or the statistical significance of such models is almost
absent in the literature. Having validated models allows us to assess
the robustness and reproducibility of our models to ensure that the
insights we get are not merely a product of chance and are therefore
generalizable.

#### Reliability

Split-half reliability assesses whether the network structure is stable
when the data is randomly divided into two halves. High reliability
means the network structure is a consistent property of the data, not
driven by a small number of idiosyncratic sequences.

``` r

network_reliability(net_tna)
#> Split-Half Reliability (1000 iterations, split = 50%)
#>   Mean Abs. Diff.     mean = 0.0361  sd = 0.0038
#>   Median Abs. Diff.   mean = 0.0269  sd = 0.0038
#>   Pearson             mean = 0.8642  sd = 0.0290
#>   Max Abs. Diff.      mean = 0.1733  sd = 0.0391
```

#### Bootstrap Analysis

Bootstrapping repeatedly draws samples from the original dataset with
replacement to estimate the model for each sample. When edges
consistently appear across the majority of the estimated models, they
are considered stable and significant. The bootstrap also provides
confidence intervals and p-values for each edge weight, offering a
quantifiable measure of uncertainty and robustness for each transition
in the network.

``` r

set.seed(42)
boot <- bootstrap_network(net_tna, iter = 100)
boot
#>   Edge                   Mean     95% CI          p
#>   -----------------------------------------------
#>   Request → Specify    0.578  [0.506, 0.661]  ** 
#>   Refine → Specify     0.417  [0.320, 0.496]  *  
#>   Command → Specify    0.294  [0.251, 0.331]  ** 
#>   Specify → Command    0.290  [0.256, 0.327]  ** 
#>   Correct → Specify    0.285  [0.234, 0.350]  *  
#>   ... and 3 more significant edges
#> 
#> Bootstrap Network  [Transition Network (relative) | directed]
#>   Iterations : 100  |  Nodes : 9
#>   Edges      : 8 significant / 81 total
#>   CI         : 95%  |  Inference: stability  |  CR [0.75, 1.25]
```

#### Centrality Stability

Centrality measures quantify the role or importance of each state in the
process. Centrality stability analysis quantifies how robust centrality
rankings are to case-dropping: the CS-coefficient is the maximum
proportion of cases that can be dropped while maintaining a correlation
of at least 0.7 with the original centrality values. A CS-coefficient
above 0.5 indicates stable rankings; below 0.25 indicates instability.

``` r

centrality_stability(net_tna, iter = 100)
#> Centrality Stability (100 iterations, threshold = 0.7)
#>   Drop proportions: 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
#> 
#>   CS-coefficients:
#>     InStrength       0.90
#>     OutStrength      0.50
#>     Betweenness      0.70
```

### Clustering

Clusters represent typical transition networks that recur across
different instances. Unlike communities, clusters involve the entire
network where groups of sequences are similarly interconnected and each
exhibits a distinct transition pattern with its own set of transition
probabilities. Identifying clusters captures the dynamics, revealing
typical behavioral strategies that learners frequently adopt.
[`build_clusters()`](https://saqr.me/Nestimate/reference/build_clusters.md)
computes pairwise dissimilarities between sequences and partitions them
into `k` groups, then builds a separate network for each cluster
(López-Pernas et al., 2025).

``` r

Cls <- build_clusters(net_tna, k = 3)
Clusters <- build_network(Cls, method = "tna")
Clusters
#> Group Networks (3 clusters via pam / hamming)
#> 
#>   Group      Nodes  Edges  Weights         N
#>   Cluster 1  9      81     [0.010, 0.626]  89 (71.2%)
#>   Cluster 2  9      80     [0.017, 0.506]  35 (28.0%)
#>   Cluster 3  9      53     [0.043, 0.667]  1 ( 0.8%)
```

#### Permutation Test for Clusters

Permutation testing is particularly important for data-driven clusters:
because clustering algorithms partition sequences to maximize
between-group separation, some degree of apparent difference is
guaranteed by construction. The permutation test provides the necessary
corrective — by randomly reassigning sequences to groups while
preserving internal sequential structure, it constructs null
distributions for edge-level differences. Only differences that exceed
this null distribution constitute evidence of genuine structural
divergence rather than algorithmic artifacts.

``` r

perm <- permutation(Clusters$`Cluster 1`, Clusters$`Cluster 2`,
                         iter = 100)
perm
#> Permutation Test:Transition Network (relative probabilities) [directed]
#>   Iterations: 100  |  Alpha: 0.05
#>   Nodes: 9  |  Edges tested: 81  |  Significant: 8
```

### Post-hoc Covariate Analysis

[`build_clusters()`](https://saqr.me/Nestimate/reference/build_clusters.md)
supports post-hoc covariate analysis: covariates do not influence the
clustering but are analyzed after the fact to characterize who ends up
in which cluster. This is the appropriate approach when the clustering
should reflect behavioral patterns alone, and the researcher then asks
whether those patterns are associated with external variables.

The `group_regulation_long` dataset contains self-regulated learning
sequences annotated with an `Achiever` covariate distinguishing high and
low achievers.

``` r

data("group_regulation_long")
net_GR <- build_network(group_regulation_long, method = "tna",
                        action = "Action", actor = "Actor",
                        time = "Time")
```

``` r

Post <- build_clusters(net_GR, k = 2, covariates = c("Achiever"))
summary(Post)
#> Sequence Clustering Summary
#>   Method:        pam 
#>   Dissimilarity: hamming 
#>   Silhouette:    0.1839 
#> 
#> Per-cluster statistics:
#>  cluster size mean_within_dist
#>        1  982         10.69340
#>        2 1018         18.59498
#> 
#> Post-hoc Covariate Analysis (does not influence cluster membership) 
#> 
#> Cluster Profiles (categorical):
#>  Cluster N    Achiever=High N(%) Achiever=Low N(%)
#>  1        982 504 (51%)          478 (49%)        
#>  2       1018 496 (49%)          522 (51%)        
#> 
#> Predictors of Membership (estimator = firth, reference: Cluster 1):
#>  Cluster Variable    OR   95% CI       p     Sig
#>  2       AchieverLow 1.11 [0.93, 1.32] 0.245    
#> 
#> Model: AIC = 2774.6 | BIC = 2785.8 | McFadden R-squared = 0.00
#> 
#> Note: Covariates are post-hoc and do not influence cluster assignments.
```

``` r

Postgr <- build_network(Post)
Postgr
#> Group Networks (2 clusters via pam / hamming)
#> 
#>   Group      Nodes  Edges  Weights         N
#>   Cluster 1  9      78     [0.001, 0.506]  982 (49.1%)
#>   Cluster 2  9      78     [0.001, 0.494]  1018 (50.9%)
```

## Part II: Psychological Network Analysis

### Theoretical Grounding

Psychological network analysis estimates the conditional dependency
structure among a set of variables. Variables (e.g., symptoms, traits,
behaviors) are represented as nodes, and edges represent partial
correlations — the association between two variables after controlling
for all others. This approach reveals which variables are directly
connected versus those whose association is mediated through other
variables (Saqr et al., 2024).

`Nestimate` supports three estimation methods for psychological
networks, all accessed through the same
[`build_network()`](https://saqr.me/Nestimate/reference/build_network.md)
interface:

- **Correlation networks** (`method = "cor"`) estimate pairwise Pearson
  correlations and produce fully connected undirected networks. While
  informative as a starting point, they do not distinguish direct from
  indirect associations.
- **Partial correlation networks** (`method = "pcor"`) control for all
  other variables, revealing only direct associations. They provide a
  more accurate picture of the dependency structure, though results can
  be noisy in small samples.
- **Regularized networks via EBICglasso** (`method = "glasso"`) apply L1
  regularization to the precision matrix, shrinking weak or unreliable
  edges to exactly zero. This is the recommended approach for
  psychological network analysis, as it balances model fit against
  complexity and produces sparse, interpretable, and replicable
  structures.

### Data

The `chatgpt_srl` dataset contains scale scores on five self-regulated
learning (SRL) constructs — Comprehension and Study Understanding (CSU),
Intrinsic Value (IV), Self-Efficacy (SE), Self-Regulation (SR), and Task
Avoidance (TA) — for 1,000 responses generated by ChatGPT to a validated
SRL questionnaire (Vogelsmeier et al., 2025).

``` r

data(chatgpt_srl)
head(chatgpt_srl)
#>        CSU       IV       SE       SR   TA
#> 1 4.538462 5.222222 4.000000 3.222222 5.00
#> 2 4.153846 4.333333 4.888889 4.333333 4.50
#> 3 4.846154 5.111111 3.777778 4.666667 4.25
#> 4 4.307692 4.222222 4.777778 3.777778 4.00
#> 5 4.307692 4.777778 4.333333 4.111111 5.75
#> 6 4.692308 3.777778 5.222222 3.888889 5.00
```

### Regularized Network (EBICglasso)

The graphical lasso applies L1 regularization to the precision matrix
(the inverse of the covariance matrix), producing a sparse network where
weak or unreliable edges are shrunk to exactly zero. The `gamma`
parameter controls sparsity through EBIC model selection — higher values
yield sparser networks.

``` r

net_glasso <- build_network(chatgpt_srl, method = "glasso",
                            params = list(gamma = 0.5))
net_glasso
#> Partial Correlation Network (EBICglasso) [undirected]
#>   Sample size: 1000
#> 
#>   Weight matrix:
#>       CSU IV SE SR TA
#>   CSU   0  0  0  0  0
#>   IV    0  0  0  0  0
#>   SE    0  0  0  0  0
#>   SR    0  0  0  0  0
#>   TA    0  0  0  0  0 
#> 
#>   Predictability (R²):
#>   CSU           0.000  
#>   IV            0.000  
#>   SE            0.000  
#>   SR            0.000  
#>   TA            0.000  
#> 
#>   Gamma: 0.50  |  Lambda: 0.0517
```

## References

Saqr, M. (2024). Temporal Network Analysis: Introduction, Methods and
Analysis with R. In M. Saqr & S. López-Pernas (Eds.), *Learning
Analytics Methods and Tutorials: A Practical Guide Using R*. Springer.
<https://lamethods.org/book1/chapters/ch17-temporal-networks/ch17-tna.html>

Saqr, M., Beck, E., & López-Pernas, S. (2024). Psychological Networks: A
Modern Approach to Analysis of Learning and Complex Learning Processes.
In M. Saqr & S. López-Pernas (Eds.), *Learning Analytics Methods and
Tutorials: A Practical Guide Using R*. Springer.
<https://lamethods.org/book1/chapters/ch19-psychological-networks/ch19-psych.html>

Saqr, M., López-Pernas, S., & Tikka, S. (2025a). Mapping Relational
Dynamics with Transition Network Analysis: A Primer and Tutorial. In M.
Saqr & S. López-Pernas (Eds.), *Advanced Learning Analytics Methods: AI,
Precision and Complexity*. Springer Nature Switzerland.
<https://lamethods.org/book2/chapters/ch15-tna/ch15-tna.html>

Saqr, M., López-Pernas, S., Törmänen, T., Kaliisa, R., Misiejuk, K., &
Tikka, S. (2025b). Transition network analysis: A novel framework for
modeling, visualizing, and identifying the temporal patterns of learners
and learning processes. *Proceedings of the 15th International Learning
Analytics and Knowledge Conference (LAK25)*. ACM.
<https://doi.org/10.1145/3706468.3706513>

López-Pernas, S., Tikka, S., & Saqr, M. (2025). Mining Patterns and
Clusters with Transition Network Analysis: A Heterogeneity Approach. In
M. Saqr & S. López-Pernas (Eds.), *Advanced Learning Analytics Methods:
AI, Precision and Complexity*. Springer Nature Switzerland.
<https://lamethods.org/book2/chapters/ch17-tna-clusters/ch17-tna-clusters.html>

Vogelsmeier, L.V.D.E., Oliveira, E., Misiejuk, K., López-Pernas, S., &
Saqr, M. (2025). Delving into the psychology of Machines: Exploring the
structure of self-regulated learning via LLM-generated survey responses.
*Computers in Human Behavior*, 173, 108769.
<https://doi.org/10.1016/j.chb.2025.108769>
