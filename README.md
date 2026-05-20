# Nestimate <img src="man/figures/logo.png" align="right" height="139"  alt="" />

> Unified network estimation, analysis, and validation for behavioral, psychological, and panel data.

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)](https://github.com/mohsaqr/Nestimate)
[![CRAN status](https://www.r-pkg.org/badges/version/Nestimate)](https://CRAN.R-project.org/package=Nestimate)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Nestimate is a comprehensive R package for estimating, validating, and comparing networks from behavioral sequence data, psychological scales, and longitudinal panel data. A single entry point — `build_network()` — dispatches to 13 built-in estimators. Every network type shares the same validation pipeline: bootstrap confidence intervals, permutation testing, split-half reliability, and centrality stability. The entire package has only 4 hard imports (ggplot2, glasso, data.table, cluster).

### Full tutorials

- [Network Estimation & Visualization](https://saqr.me/Nestimate/articles/cograph-tutorial-nestimate.html)
- [Higher-Order & Simplicial Complexes](https://saqr.me/Nestimate/articles/cograph-tutorial-simplicial.html)
- [Model Assessment](https://saqr.me/Nestimate/articles/tutorial_model_assessment.html)

### Quick guides

- [Transition Networks](https://saqr.me/Nestimate/articles/transition-networks.html)
- [Sequence Plots & Comparison](https://saqr.me/Nestimate/articles/sequence-plots.html)
- [Clustering & Multi-Level Analysis](https://saqr.me/Nestimate/articles/clustering.html)
- [Markov Stability](https://saqr.me/Nestimate/articles/markov-stability.html)
- [Sequence Pattern Comparison](https://saqr.me/Nestimate/articles/sequence-comparison.html)


## Installation

```r
# From CRAN
install.packages("Nestimate")

# Development version
devtools::install_github("mohsaqr/Nestimate")
```

## What Nestimate Covers

| Area | Key Functions |
|------|--------------|
| Dynamic / Transition Networks | `build_network()`, `wtna()`, `cooccurrence()` |
| Psychological Networks | `build_network(method = "glasso/pcor/cor/ising/mgm")` |
| Multilevel VAR | `build_mlvar()` |
| Idiographic Networks | `build_gimme()` |
| Cluster & Group Networks | `build_clusters()`, `build_mmm()`, `build_mcml()` |
| Higher-Order Networks | `build_hon()`, `build_honem()`, `build_hypa()`, `build_mogen()` |
| Topological Analysis | `build_simplicial()`, `persistent_homology()`, `q_analysis()` |
| Sequence Visualization | `sequence_plot()`, `distribution_plot()` |
| Sequence Pattern Comparison | `sequence_compare()` |
| Association Mining | `association_rules()` |
| Link Prediction | `predict_links()`, `evaluate_links()` |
| Markov Chain Analysis | `markov_stability()`, `passage_time()` |
| Statistical Validation | `bootstrap_network()`, `permutation()`, `nct()`, `network_reliability()`, `centrality_stability()` |

---

## Dynamic Networks

All dynamic network methods use `build_network()`. Pass an event log with `action`, `actor`, and `time` columns — no preprocessing needed.

### Estimation Methods

| Method | Aliases | Description |
|--------|---------|-------------|
| `"relative"` | `"tna"`, `"transition"` | Transition probabilities (directed) |
| `"frequency"` | `"ftna"`, `"counts"` | Raw transition counts (directed) |
| `"attention"` | `"atna"` | Decay-weighted transitions emphasising recent events (directed) |
| `"co_occurrence"` | `"cna"` | Co-occurrence from sequential data (undirected) |

```r
library(Nestimate)
data(human_long)

net <- build_network(human_long, method = "tna",
                     action = "action", actor = "session_id", time = "time")

# Per-group networks in one call
group_nets <- build_network(human_long, method = "tna",
                            action = "action", actor = "session_id",
                            time = "time", group = "phase")
```

### Window-Based TNA

`wtna()` builds networks from binary (one-hot) data using temporal windows — directed transitions between windows, undirected co-occurrence within windows, or a mixed network combining both.

```r
data(learning_activities)
net_wtna  <- wtna(learning_activities, actor = "student",
                  method = "transition", type = "relative")
net_mixed <- wtna(learning_activities, actor = "student",
                  method = "both", type = "relative")
```

### Co-occurrence Networks

`cooccurrence()` builds undirected co-occurrence networks from 6 input formats (delimited fields, long/bipartite, binary matrix, wide sequence, lists) with 8 similarity methods (Jaccard, cosine, association strength, Dice, and more).

```r
# From a long-format data frame
net_co <- cooccurrence(human_long, field = "action", by = "session_id",
                       similarity = "jaccard", threshold = 0.1)
```

---

## Psychological Networks

| Method | Description |
|--------|-------------|
| `"cor"` | Pearson correlations |
| `"pcor"` | Partial correlations (precision matrix inversion) |
| `"glasso"` | EBICglasso — sparse regularised partial correlations |
| `"ising"` | L1-regularised logistic regression for binary items |
| `"mgm"` | Mixed Graphical Model — continuous + categorical variables together |

All implemented from scratch with no dependency on igraph, qgraph, or bootnet.

```r
data(srl_strategies)
net_gl  <- build_network(srl_strategies, method = "glasso")
net_mgm <- build_network(mixed_data, method = "mgm")   # scales + demographics

predictability(net_gl)   # R-squared per node from network structure
```

---

## Multilevel VAR

`build_mlvar()` estimates three networks simultaneously from ESM/EMA diary data — the three pillars of mlVAR analysis in a single function call:

- **Temporal** (directed) — lagged fixed effects: what predicts what across time
- **Contemporaneous** (undirected) — partial correlations of within-person residuals: what moves together right now
- **Between-subjects** (undirected) — partial correlations of person means: who differs from whom

Machine-precision equivalence to `mlVAR::mlVAR()` validated across 25 real ESM datasets, runs 1.45× faster.

```r
data(chatgpt_srl)
fit <- build_mlvar(chatgpt_srl, vars = c("planning", "monitoring", "evaluation"),
                   id = "id", day = "day", beep = "beep")

fit$temporal          # directed network of lagged effects
fit$contemporaneous   # undirected within-person partial correlations
fit$between           # undirected between-persons partial correlations
coefs(fit)            # tidy data.frame: beta, SE, t, p, CI for every edge
```

---

## Idiographic Networks

`build_gimme()` estimates a separate network for each person using the Group Iterative Mean Estimation (GIMME) algorithm, then aggregates to a group-level picture. Use this when between-person heterogeneity matters and a single group network would average over meaningfully different individuals.

```r
fit_g <- build_gimme(panel_data, vars = c("x1", "x2", "x3"), id = "id")
fit_g$group_network       # aggregated group-level paths
fit_g$individual_networks # one network per person
```

---

## Cluster & Group Networks

### Sequence Clustering

`build_clusters()` partitions sequences into `k` groups using pairwise distance matrices. Supports 9 distance metrics and 8 clustering algorithms. Both `build_clusters()` and `build_mmm()` results pass directly to `build_network()`.

```r
clust <- build_clusters(net, k = 3, dissimilarity = "hamming", method = "ward.D2")
plot(clust, type = "silhouette")
cluster_nets <- build_network(clust, method = "tna")
```

### Mixed Markov Models

`build_mmm()` discovers latent subgroups of sequences that share similar transition dynamics via EM — without pre-labelling groups. BIC/AIC/ICL model selection via `compare_mmm()`.

```r
mmm   <- build_mmm(net, k = 3)
compare_mmm(net, k = 2:6)   # model selection plot + table
mmm_nets <- build_network(mmm, method = "tna")
```

### MCML

`build_mcml()` decomposes a network into macro (between-cluster) and micro (within-cluster) layers when nodes belong to known groups.

```r
clusters <- list(Metacognitive = c("Planning", "Monitoring"),
                 Cognitive     = c("Elaboration", "Organisation"))
mcml <- cluster_summary(net, clusters)
mcml$macro$weights
```

---

## Higher-Order Networks

Capture dependencies beyond first-order transitions:

| Function | What it finds |
|----------|--------------|
| `build_hon()` | Variable-length memory paths |
| `build_honem()` | Higher-order network embedding |
| `build_hypa()` | Statistically anomalous paths (over/under-represented) |
| `build_mogen()` | Optimal Markov order per node |

```r
hon  <- build_hon(net, max_order = 2)
pathways(hon)                      # arrow-notation path strings
hypa <- build_hypa(net)
hypa$over                          # over-represented paths with p-values
```

---

## Topological Analysis

Go beyond edges — find cliques, holes, and high-order connectivity using tools from algebraic topology.

```r
sc <- build_simplicial(net, method = "clique")
betti_numbers(sc)          # connected components, cycles, voids
euler_characteristic(sc)
ph <- persistent_homology(net)  # track topology across thresholds
plot(ph)
qa <- q_analysis(sc)       # Atkin's Q-connectivity structure vectors
```

---

## Sequence Visualization

Visualize raw sequence data as index plots or state distribution charts — before or after clustering.

```r
# Sequence index plot: one row per person, coloured by state
sequence_plot(net)

# After clustering: faceted by cluster
clust <- build_clusters(net, k = 3)
sequence_plot(clust, type = "index")

# State distribution over time
distribution_plot(net, type = "area")
distribution_plot(clust, type = "bar")
```

---

## Sequence Pattern Comparison

`sequence_compare()` extracts all k-gram patterns from grouped sequences, counts per-group frequencies, and tests statistical differences via permutation — answering the question "do these groups actually behave differently, and where?"

```r
data(human_long)
net <- build_network(human_long, method = "tna",
                     action = "action", actor = "session_id",
                     time = "time", group = "phase")

res <- sequence_compare(net, sub = 2:4, test = "chisq", adjust = "fdr")
res$patterns                        # per-pattern frequencies + p-values
plot(res)                           # back-to-back pyramid chart
plot(res, style = "heatmap")        # heatmap for many patterns
```

---

## Association Rule Mining

`association_rules()` mines "if A then B" patterns from sequences or binary matrices using the Apriori algorithm. Returns support, confidence, lift, and conviction for every rule above a threshold.

```r
rules <- association_rules(net, min_support = 0.05, min_confidence = 0.6)
rules$rules                       # tidy data.frame, sorted by lift
pathways(rules)                   # rules as arrow-notation strings

# From a raw binary matrix
rules2 <- association_rules(binary_mat, min_support = 0.1)
```

---

## Link Prediction

`predict_links()` scores all unobserved node pairs using structural similarity, identifying which missing connections are most likely to exist. `evaluate_links()` computes AUC, precision, and recall against held-out edges.

```r
preds <- predict_links(net)       # common neighbours, Adamic-Adar, Katz, ...
head(preds$scores)                # sorted by predicted score

# Evaluate against known missing edges
eval  <- evaluate_links(net, held_out = test_edges)
eval$auc
```

---

## Markov Chain Analysis

`markov_stability()` measures how stable a network partition is under random-walk dynamics at different time scales — a resolution-free way to find communities. `passage_time()` computes expected first-passage and return times between states.

```r
stab <- markov_stability(net, times = seq(0.1, 10, 0.1))
plot(stab)                # stability vs time-scale curve

pt <- passage_time(net)
pt$first_passage          # expected steps to reach state j from state i
pt$return_time            # expected steps to return to the same state
```

---

## Statistical Validation

Every network type shares the same validation pipeline.

```r
# Bootstrap confidence intervals
boot <- bootstrap_network(net, iter = 1000)
summary(boot)

# Permutation test: are two networks different?
perm <- permutation(group_nets$`Cluster 1`, group_nets$`Cluster 2`)

# Network Comparison Test (NCT): formal test of structure + global strength
nct_res <- nct(data1, data2, iter = 500)
print(nct_res)            # M-statistic, S-statistic, per-edge p-values

# Split-half reliability
network_reliability(net)

# Centrality stability (CS-coefficient)
centrality_stability(net)

# Glasso-specific bootstrap (edge inclusion + centrality CIs)
boot_gl <- boot_glasso(net_pna, iter = 1000)
```

| Function | Purpose |
|----------|---------|
| `bootstrap_network()` | Bootstrap CIs and p-values for each edge |
| `permutation()` | Edge-level comparison between two networks |
| `nct()` | Formal Network Comparison Test (global strength + structure) |
| `network_reliability()` | Split-half reliability of edge weights |
| `centrality_stability()` | CS-coefficient via case-dropping |
| `boot_glasso()` | Edge inclusion, centrality CIs, difference tests for glasso networks |

---

## Bundled Datasets

| Dataset | Description |
|---------|-------------|
| `human_long` | 10,796 human actions across 429 human-AI coding sessions |
| `ai_long` | Matched AI actions from the same 429 sessions |
| `srl_strategies` | SRL strategy frequencies — 250 students, 9 strategies |
| `chatgpt_srl` | ChatGPT-generated SRL scale scores for psychological networks |
| `learning_activities` | Binary learning activity indicators — 200 students × 30 timepoints |
| `group_regulation_long` | Group regulation sequences with covariates |
| `trajectories` | 138-student engagement trajectory matrix |

---

## Documentation

- [Transition Networks](https://saqr.me/Nestimate/articles/transition-networks.html)
- [Clustering & Multilevel](https://saqr.me/Nestimate/articles/clustering.html)
- [Full Reference](https://saqr.me/Nestimate/reference/)

## Citation

If you use Nestimate in your research, please cite:

> Saqr, M., Lopez-Pernas, S., Tormanen, T., Kaliisa, R., Misiejuk, K., & Tikka, S. (2025). Transition Network Analysis: A Novel Framework for Modeling, Visualizing, and Identifying the Temporal Patterns of Learners and Learning. *Proceedings of the 15th Learning Analytics and Knowledge Conference*. doi: [10.1145/3706468.3706513](https://doi.org/10.1145/3706468.3706513)

> Saqr, M., Beck, E., & Lopez-Pernas, S. (2024). Psychological Networks. In M. Saqr & S. Lopez-Pernas (Eds.), *Learning Analytics Methods and Tutorials* (pp. 513–546). Springer. doi: [10.1007/978-3-031-54464-4_19](https://doi.org/10.1007/978-3-031-54464-4_19)

## License

MIT
