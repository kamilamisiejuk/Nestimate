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
head(human_long)
#>   message_id   project   session_id  timestamp session_date      code
#> 1       3439 Project_7 0086cabebd15 1772661600   2026-03-05   Specify
#> 2       3439 Project_7 0086cabebd15 1772661600   2026-03-05   Command
#> 3       3439 Project_7 0086cabebd15 1772661600   2026-03-05   Specify
#> 4       3440 Project_7 0086cabebd15 1772661600   2026-03-05 Interrupt
#> 5       3442 Project_7 0086cabebd15 1772661600   2026-03-05    Verify
#> 6       3444 Project_7 0086cabebd15 1772661600   2026-03-05   Specify
#>         cluster code_order order_in_session
#> 1     Directive          1                1
#> 2     Directive          2                2
#> 3     Directive          3                3
#> 4 Metacognitive          1                4
#> 5    Evaluative          1                7
#> 6     Directive          1               10
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

net_tna <- build_network(human_long, method = "tna",
                         action = "code", actor = "session_id",
                         time = "timestamp")
print(net_tna)
#> Transition Network (relative probabilities) [directed]
#>   Weights: [0.018, 0.620]  |  mean: 0.111
#> 
#>   Weight matrix:
#>             Command Correct Frustrate Inquire Interrupt Refine Request Specify
#>   Command     0.235   0.088     0.055   0.066     0.035  0.038   0.155   0.280
#>   Correct     0.093   0.091     0.138   0.057     0.054  0.112   0.120   0.278
#>   Frustrate   0.103   0.115     0.176   0.075     0.047  0.171   0.111   0.125
#>   Inquire     0.196   0.126     0.098   0.188     0.094  0.062   0.084   0.106
#>   Interrupt   0.259   0.081     0.094   0.123     0.102  0.080   0.069   0.122
#>   Refine      0.058   0.075     0.072   0.047     0.042  0.086   0.146   0.457
#>   Request     0.096   0.019     0.044   0.067     0.051  0.033   0.039   0.620
#>   Specify     0.263   0.058     0.081   0.072     0.167  0.074   0.070   0.172
#>   Verify      0.224   0.079     0.164   0.118     0.043  0.097   0.116   0.083
#>             Verify
#>   Command    0.048
#>   Correct    0.056
#>   Frustrate  0.077
#>   Inquire    0.044
#>   Interrupt  0.069
#>   Refine     0.018
#>   Request    0.032
#>   Specify    0.043
#>   Verify     0.077 
#> 
#>   Initial probabilities:
#>   Specify       0.679  ████████████████████████████████████████
#>   Command       0.175  ██████████
#>   Request       0.055  ███
#>   Interrupt     0.030  ██
#>   Correct       0.019  █
#>   Refine        0.019  █
#>   Frustrate     0.013  █
#>   Inquire       0.008  
#>   Verify        0.002
```

#### Frequency Network (FTNA)

The frequency method preserves raw transition counts rather than
normalizing to conditional probabilities. This is useful when absolute
frequencies matter — a transition occurring 500 times from a common
state may be more practically important than one occurring 5 times from
a rare state, even if the latter has a higher conditional probability.

``` r

net_ftna <- build_network(human_long, method = "ftna",
                          action = "code", actor = "session_id",
                          time = "timestamp")
print(net_ftna)
#> Transition Network (frequency counts) [directed]
#>   Weights: [14.000, 732.000]  |  mean: 126.790
#> 
#>   Weight matrix:
#>             Command Correct Frustrate Inquire Interrupt Refine Request Specify
#>   Command       459     171       107     129        68     75     303     548
#>   Correct        73      71       108      45        42     88      94     218
#>   Frustrate      95     106       162      69        43    157     102     115
#>   Inquire       160     103        80     154        77     51      69      87
#>   Interrupt     191      60        69      91        75     59      51      90
#>   Refine         45      59        56      37        33     67     114     357
#>   Request        97      19        45      68        52     33      40     629
#>   Specify       732     161       226     199       463    205     194     479
#>   Verify        108      38        79      57        21     47      56      40
#>             Verify
#>   Command       94
#>   Correct       44
#>   Frustrate     71
#>   Inquire       36
#>   Interrupt     51
#>   Refine        14
#>   Request       32
#>   Specify      120
#>   Verify        37 
#> 
#>   Initial probabilities:
#>   Specify       0.679  ████████████████████████████████████████
#>   Command       0.175  ██████████
#>   Request       0.055  ███
#>   Interrupt     0.030  ██
#>   Correct       0.019  █
#>   Refine        0.019  █
#>   Frustrate     0.013  █
#>   Inquire       0.008  
#>   Verify        0.002
```

#### Attention Network (ATNA)

The attention method applies temporal decay weighting, giving more
importance to recent transitions within each sequence. The `lambda`
parameter controls the decay rate: higher values produce faster decay.
This captures the idea that later events in a process may be more
indicative of the underlying dynamics than early ones.

``` r

net_atna <- build_network(human_long, method = "atna",
                          action = "code", actor = "session_id",
                          time = "timestamp")
print(net_atna)
#> Attention Network (decay-weighted transitions) [directed]
#>   Weights: [11.018, 357.662]  |  mean: 71.690
#> 
#>   Weight matrix:
#>             Command Correct Frustrate Inquire Interrupt  Refine Request Specify
#>   Command   254.753  91.496    66.797  80.402    65.936  47.423 145.985 297.131
#>   Correct    51.382  44.627    55.538  29.832    27.286  46.929  46.160 119.340
#>   Frustrate  57.452  55.341    92.756  40.085    28.624  74.892  54.606  87.031
#>   Inquire    85.134  56.787    43.046  81.061    38.782  29.220  40.091  69.990
#>   Interrupt  94.228  31.160    36.365  46.634    48.698  30.005  28.457  67.398
#>   Refine     36.306  37.008    37.761  24.197    21.236  51.087  58.478 167.226
#>   Request    63.272  18.659    31.998  41.056    51.503  30.285  33.277 272.476
#>   Specify   357.662  96.031   134.432 115.671   206.063 116.889 141.240 335.302
#>   Verify     53.941  20.958    40.429  30.820    11.349  25.768  29.387  39.412
#>             Verify
#>   Command   53.492
#>   Correct   22.847
#>   Frustrate 32.765
#>   Inquire   17.833
#>   Interrupt 27.393
#>   Refine    11.018
#>   Request   21.781
#>   Specify   76.312
#>   Verify    23.444 
#> 
#>   Initial probabilities:
#>   Specify       0.679  ████████████████████████████████████████
#>   Command       0.175  ██████████
#>   Request       0.055  ███
#>   Interrupt     0.030  ██
#>   Correct       0.019  █
#>   Refine        0.019  █
#>   Frustrate     0.013  █
#>   Inquire       0.008  
#>   Verify        0.002
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
#>   Mean Abs. Diff.     mean = 0.0169  sd = 0.0017
#>   Median Abs. Diff.   mean = 0.0129  sd = 0.0018
#>   Pearson             mean = 0.9706  sd = 0.0062
#>   Max Abs. Diff.      mean = 0.0748  sd = 0.0181
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
#>   Request → Specify    0.619  [0.584, 0.655]  ** 
#>   Refine → Specify     0.455  [0.413, 0.502]  ** 
#>   Command → Specify    0.281  [0.256, 0.301]  ** 
#>   Correct → Specify    0.280  [0.249, 0.309]  ** 
#>   Specify → Command    0.263  [0.246, 0.279]  ** 
#>   ... and 40 more significant edges
#> 
#> Bootstrap Network  [Transition Network (relative) | directed]
#>   Iterations : 100  |  Nodes : 9
#>   Edges      : 45 significant / 81 total
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
#>     OutStrength      0.80
#>     Betweenness      0.90
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
#>   Cluster 1  9      81     [0.006, 0.538]  145 (27.6%)
#>   Cluster 2  9      81     [0.003, 0.694]  299 (56.8%)
#>   Cluster 3  9      81     [0.018, 0.588]  82 (15.6%)
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
#>   Nodes: 9  |  Edges tested: 81  |  Significant: 9
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
