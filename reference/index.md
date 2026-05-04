# Package index

## Network Estimation

Core functions for building networks from data

- [`build_network()`](https://mohsaqr.github.io/Nestimate/reference/build_network.md)
  : Build a Network
- [`estimate_network()`](https://mohsaqr.github.io/Nestimate/reference/estimate_network.md)
  : Estimate a Network (Deprecated)
- [`register_estimator()`](https://mohsaqr.github.io/Nestimate/reference/register_estimator.md)
  : Register a Network Estimator
- [`get_estimator()`](https://mohsaqr.github.io/Nestimate/reference/get_estimator.md)
  : Retrieve a Registered Estimator
- [`list_estimators()`](https://mohsaqr.github.io/Nestimate/reference/list_estimators.md)
  : List All Registered Estimators
- [`remove_estimator()`](https://mohsaqr.github.io/Nestimate/reference/remove_estimator.md)
  : Remove a Registered Estimator
- [`build_tna()`](https://mohsaqr.github.io/Nestimate/reference/build_tna.md)
  : Build a Transition Network (TNA)
- [`build_atna()`](https://mohsaqr.github.io/Nestimate/reference/build_atna.md)
  : Build an Attention-Weighted Transition Network (ATNA)
- [`build_ftna()`](https://mohsaqr.github.io/Nestimate/reference/build_ftna.md)
  : Build a Frequency Transition Network (FTNA)
- [`build_cna()`](https://mohsaqr.github.io/Nestimate/reference/build_cna.md)
  : Build a Co-occurrence Network (CNA)
- [`build_cor()`](https://mohsaqr.github.io/Nestimate/reference/build_cor.md)
  : Build a Correlation Network
- [`build_pcor()`](https://mohsaqr.github.io/Nestimate/reference/build_pcor.md)
  : Build a Partial Correlation Network
- [`build_glasso()`](https://mohsaqr.github.io/Nestimate/reference/build_glasso.md)
  : Build a Graphical Lasso Network (EBICglasso)
- [`build_ising()`](https://mohsaqr.github.io/Nestimate/reference/build_ising.md)
  : Build an Ising Network
- [`wtna()`](https://mohsaqr.github.io/Nestimate/reference/wtna.md) :
  Window-based Transition Network Analysis
- [`cooccurrence()`](https://mohsaqr.github.io/Nestimate/reference/cooccurrence.md)
  : Build a Co-occurrence Network
- [`build_mlvar()`](https://mohsaqr.github.io/Nestimate/reference/build_mlvar.md)
  : Build a Multilevel Vector Autoregression (mlVAR) network
- [`build_gimme()`](https://mohsaqr.github.io/Nestimate/reference/build_gimme.md)
  : GIMME: Group Iterative Multiple Model Estimation

## Higher-Order Networks

Methods for capturing higher-order dependencies

- [`build_hon()`](https://mohsaqr.github.io/Nestimate/reference/build_hon.md)
  : Build a Higher-Order Network (HON)
- [`build_honem()`](https://mohsaqr.github.io/Nestimate/reference/build_honem.md)
  : Build HONEM Embeddings for Higher-Order Networks
- [`build_hypa()`](https://mohsaqr.github.io/Nestimate/reference/build_hypa.md)
  : Detect Path Anomalies via HYPA
- [`build_mogen()`](https://mohsaqr.github.io/Nestimate/reference/build_mogen.md)
  : Build Multi-Order Generative Model (MOGen)
- [`pathways()`](https://mohsaqr.github.io/Nestimate/reference/pathways.md)
  : Extract Pathways from Higher-Order Network Objects
- [`mogen_transitions()`](https://mohsaqr.github.io/Nestimate/reference/mogen_transitions.md)
  : Extract Transition Table from a MOGen Model
- [`path_counts()`](https://mohsaqr.github.io/Nestimate/reference/path_counts.md)
  : Count Path Frequencies in Trajectory Data
- [`bipartite_groups()`](https://mohsaqr.github.io/Nestimate/reference/bipartite_groups.md)
  : Hypergraph from bipartite group / event data
- [`clique_expansion()`](https://mohsaqr.github.io/Nestimate/reference/clique_expansion.md)
  : Clique expansion of a hypergraph
- [`hypergraph_centrality()`](https://mohsaqr.github.io/Nestimate/reference/hypergraph_centrality.md)
  : Hypergraph eigenvector centralities

## Markov Analysis

Order, structure, entropy, and stability of Markov chains

- [`chain_structure()`](https://mohsaqr.github.io/Nestimate/reference/chain_structure.md)
  : Qualitative structure of a discrete-time Markov chain.
- [`markov_order_test()`](https://mohsaqr.github.io/Nestimate/reference/markov_order_test.md)
  : Test the Markov order of a sequential process
- [`markov_stability()`](https://mohsaqr.github.io/Nestimate/reference/markov_stability.md)
  [`plot(`*`<net_markov_stability>`*`)`](https://mohsaqr.github.io/Nestimate/reference/markov_stability.md)
  : Markov Stability Analysis
- [`passage_time()`](https://mohsaqr.github.io/Nestimate/reference/passage_time.md)
  [`summary(`*`<net_mpt>`*`)`](https://mohsaqr.github.io/Nestimate/reference/passage_time.md)
  [`plot(`*`<net_mpt>`*`)`](https://mohsaqr.github.io/Nestimate/reference/passage_time.md)
  : Mean First Passage Times
- [`path_dependence()`](https://mohsaqr.github.io/Nestimate/reference/path_dependence.md)
  : Per-Context Path Dependence at Order k
- [`transition_entropy()`](https://mohsaqr.github.io/Nestimate/reference/transition_entropy.md)
  : Transition Entropy of a Markov Chain

## Bootstrap & Inference

Statistical inference for network estimation

- [`bootstrap_network()`](https://mohsaqr.github.io/Nestimate/reference/bootstrap_network.md)
  : Bootstrap a Network Estimate
- [`boot_glasso()`](https://mohsaqr.github.io/Nestimate/reference/boot_glasso.md)
  : Bootstrap for Regularized Partial Correlation Networks
- [`permutation()`](https://mohsaqr.github.io/Nestimate/reference/permutation.md)
  : Permutation Test for Network Comparison
- [`nct()`](https://mohsaqr.github.io/Nestimate/reference/nct.md) :
  Network Comparison Test

## Reliability & Stability

Assess reliability and stability of network estimates

- [`network_reliability()`](https://mohsaqr.github.io/Nestimate/reference/network_reliability.md)
  : Split-Half Reliability for Network Estimates
- [`centrality_stability()`](https://mohsaqr.github.io/Nestimate/reference/centrality_stability.md)
  : Centrality Stability Coefficient (CS-coefficient)

## Clustering & Grouping

Cluster-based and multilevel network analysis

- [`build_clusters()`](https://mohsaqr.github.io/Nestimate/reference/build_clusters.md)
  : Cluster Sequences by Dissimilarity
- [`cluster_choice()`](https://mohsaqr.github.io/Nestimate/reference/cluster_choice.md)
  : Cluster Choice – sweep k, dissimilarity and method
- [`cluster_diagnostics()`](https://mohsaqr.github.io/Nestimate/reference/cluster_diagnostics.md)
  [`as.data.frame(`*`<net_cluster_diagnostics>`*`)`](https://mohsaqr.github.io/Nestimate/reference/cluster_diagnostics.md)
  : Cluster Diagnostics
- [`cluster_summary()`](https://mohsaqr.github.io/Nestimate/reference/cluster_summary.md)
  : Cluster Summary Statistics
- [`cluster_mmm()`](https://mohsaqr.github.io/Nestimate/reference/cluster_mmm.md)
  : Cluster sequences using Mixed Markov Models
- [`cluster_network()`](https://mohsaqr.github.io/Nestimate/reference/cluster_network.md)
  : Cluster data and build per-cluster networks in one step
- [`build_mcml()`](https://mohsaqr.github.io/Nestimate/reference/build_mcml.md)
  : Build MCML from Raw Transition Data
- [`build_mmm()`](https://mohsaqr.github.io/Nestimate/reference/build_mmm.md)
  : Fit a Mixed Markov Model
- [`compare_mmm()`](https://mohsaqr.github.io/Nestimate/reference/compare_mmm.md)
  : Compare MMM fits across different k

## Simplicial Complex Analysis

Topological analysis of networks

- [`build_simplicial()`](https://mohsaqr.github.io/Nestimate/reference/build_simplicial.md)
  : Build a Simplicial Complex
- [`persistent_homology()`](https://mohsaqr.github.io/Nestimate/reference/persistent_homology.md)
  : Persistent Homology
- [`q_analysis()`](https://mohsaqr.github.io/Nestimate/reference/q_analysis.md)
  : Q-Analysis
- [`betti_numbers()`](https://mohsaqr.github.io/Nestimate/reference/betti_numbers.md)
  : Betti Numbers
- [`euler_characteristic()`](https://mohsaqr.github.io/Nestimate/reference/euler_characteristic.md)
  : Euler Characteristic
- [`simplicial_degree()`](https://mohsaqr.github.io/Nestimate/reference/simplicial_degree.md)
  : Simplicial Degree
- [`verify_simplicial()`](https://mohsaqr.github.io/Nestimate/reference/verify_simplicial.md)
  : Verify Simplicial Complex Against igraph

## Data Preparation

Convert and prepare data for network estimation

- [`prepare()`](https://mohsaqr.github.io/Nestimate/reference/prepare.md)
  : Prepare Event Log Data for Network Estimation
- [`prepare_for_tna()`](https://mohsaqr.github.io/Nestimate/reference/prepare_for_tna.md)
  : Prepare Data for TNA Analysis
- [`action_to_onehot()`](https://mohsaqr.github.io/Nestimate/reference/action_to_onehot.md)
  : Convert Action Column to One-Hot Encoding
- [`prepare_onehot()`](https://mohsaqr.github.io/Nestimate/reference/prepare_onehot.md)
  : Import One-Hot Encoded Data into Sequence Format
- [`wide_to_long()`](https://mohsaqr.github.io/Nestimate/reference/wide_to_long.md)
  : Convert Wide Sequences to Long Format
- [`long_to_wide()`](https://mohsaqr.github.io/Nestimate/reference/long_to_wide.md)
  : Convert Long Format to Wide Sequences
- [`convert_sequence_format()`](https://mohsaqr.github.io/Nestimate/reference/convert_sequence_format.md)
  : Convert Sequence Data to Different Formats
- [`actor_endpoints()`](https://mohsaqr.github.io/Nestimate/reference/actor_endpoints.md)
  : Tidy per-actor endpoint summary of a wide-format sequence dataset.
- [`mark_first_state()`](https://mohsaqr.github.io/Nestimate/reference/mark_first_state.md)
  : Mark leading-NA cells with an explicit state label.
- [`mark_terminal_state()`](https://mohsaqr.github.io/Nestimate/reference/mark_terminal_state.md)
  : Mark terminal-NA cells with an explicit state label.

## Utilities

Helper functions and extractors

- [`predictability()`](https://mohsaqr.github.io/Nestimate/reference/predictability.md)
  : Compute Node Predictability
- [`frequencies()`](https://mohsaqr.github.io/Nestimate/reference/frequencies.md)
  : Sequence Data Conversion Functions
- [`state_frequencies()`](https://mohsaqr.github.io/Nestimate/reference/state_frequencies.md)
  : Compute State Frequencies from Trajectory Data
- [`net_aggregate_weights()`](https://mohsaqr.github.io/Nestimate/reference/net_aggregate_weights.md)
  : Aggregate Edge Weights
- [`net_centrality()`](https://mohsaqr.github.io/Nestimate/reference/net_centrality.md)
  : Compute Centrality Measures for a Network
- [`coefs()`](https://mohsaqr.github.io/Nestimate/reference/coefs.md) :
  Tidy coefficients from a fitted mlvar model
- [`as_tna()`](https://mohsaqr.github.io/Nestimate/reference/as_tna.md)
  : Convert cluster_summary to tna Objects
- [`extract_edges()`](https://mohsaqr.github.io/Nestimate/reference/extract_edges.md)
  : Extract Edge List with Weights
- [`extract_initial_probs()`](https://mohsaqr.github.io/Nestimate/reference/extract_initial_probs.md)
  : Extract Initial Probabilities from Model
- [`extract_transition_matrix()`](https://mohsaqr.github.io/Nestimate/reference/extract_transition_matrix.md)
  : Extract Transition Matrix from Model

## Sequence Analysis

Sequence visualization, pattern comparison, and association mining

- [`sequence_plot()`](https://mohsaqr.github.io/Nestimate/reference/sequence_plot.md)
  : Sequence Plot (heatmap, index, or distribution)
- [`distribution_plot()`](https://mohsaqr.github.io/Nestimate/reference/distribution_plot.md)
  : State Distribution Plot Over Time
- [`plot_state_frequencies()`](https://mohsaqr.github.io/Nestimate/reference/plot_state_frequencies.md)
  : Plot State Frequency Distributions
- [`state_distribution()`](https://mohsaqr.github.io/Nestimate/reference/state_distribution.md)
  : Per-Class State Distribution as a Tidy Data Frame
- [`plot_mosaic()`](https://mohsaqr.github.io/Nestimate/reference/plot_mosaic.md)
  : Draw a Marimekko / Mosaic Plot from a Tidy Data Frame
- [`mosaic_plot()`](https://mohsaqr.github.io/Nestimate/reference/mosaic_plot.md)
  : Mosaic Plot of a Network's Transition or Co-occurrence Counts
- [`sequence_compare()`](https://mohsaqr.github.io/Nestimate/reference/sequence_compare.md)
  : Compare Subsequence Patterns Between Groups
- [`association_rules()`](https://mohsaqr.github.io/Nestimate/reference/association_rules.md)
  : Discover Association Rules from Sequential or Transaction Data

## Link Prediction

Predict and evaluate missing connections

- [`predict_links()`](https://mohsaqr.github.io/Nestimate/reference/predict_links.md)
  : Predict Missing or Future Links in a Network
- [`evaluate_links()`](https://mohsaqr.github.io/Nestimate/reference/evaluate_links.md)
  : Evaluate Link Predictions Against Known Edges

## Data

Example datasets

- [`human_long`](https://mohsaqr.github.io/Nestimate/reference/long-data.md)
  [`ai_long`](https://mohsaqr.github.io/Nestimate/reference/long-data.md)
  : Human-AI Vibe Coding Interaction Data (Long Format)
- [`srl_strategies`](https://mohsaqr.github.io/Nestimate/reference/srl_strategies.md)
  : Self-Regulated Learning Strategy Frequencies
- [`learning_activities`](https://mohsaqr.github.io/Nestimate/reference/learning_activities.md)
  : Online Learning Activity Indicators
- [`group_regulation_long`](https://mohsaqr.github.io/Nestimate/reference/group_regulation_long.md)
  : Group Regulation in Collaborative Learning (Long Format)
- [`chatgpt_srl`](https://mohsaqr.github.io/Nestimate/reference/chatgpt_srl.md)
  : ChatGPT Self-Regulated Learning Scale Scores
- [`trajectories`](https://mohsaqr.github.io/Nestimate/reference/trajectories.md)
  : Student Engagement Trajectories

## S3 Methods

Print, summary, and plot methods for package classes

- [`build_hypergraph()`](https://mohsaqr.github.io/Nestimate/reference/build_hypergraph.md)
  [`print(`*`<net_hypergraph>`*`)`](https://mohsaqr.github.io/Nestimate/reference/build_hypergraph.md)
  [`summary(`*`<net_hypergraph>`*`)`](https://mohsaqr.github.io/Nestimate/reference/build_hypergraph.md)
  : Higher-order hypergraph from a network's clique structure

- [`casedrop_reliability()`](https://mohsaqr.github.io/Nestimate/reference/casedrop_reliability.md)
  [`print(`*`<net_casedrop_reliability>`*`)`](https://mohsaqr.github.io/Nestimate/reference/casedrop_reliability.md)
  [`summary(`*`<net_casedrop_reliability>`*`)`](https://mohsaqr.github.io/Nestimate/reference/casedrop_reliability.md)
  [`print(`*`<net_casedrop_reliability_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/casedrop_reliability.md)
  [`summary(`*`<net_casedrop_reliability_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/casedrop_reliability.md)
  [`print(`*`<summary.net_casedrop_reliability_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/casedrop_reliability.md)
  [`plot(`*`<net_casedrop_reliability>`*`)`](https://mohsaqr.github.io/Nestimate/reference/casedrop_reliability.md)
  [`plot(`*`<net_casedrop_reliability_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/casedrop_reliability.md)
  : Edge-weight Case-dropping Stability

- [`hypergraph_measures()`](https://mohsaqr.github.io/Nestimate/reference/hypergraph_measures.md)
  [`print(`*`<hypergraph_measures>`*`)`](https://mohsaqr.github.io/Nestimate/reference/hypergraph_measures.md)
  : Structural measures for a hypergraph

- [`print(`*`<boot_glasso>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.boot_glasso.md)
  : Print Method for boot_glasso

- [`print(`*`<chain_structure>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.chain_structure.md)
  :

  Print method for `chain_structure`.

- [`print(`*`<chain_structure_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.chain_structure_group.md)
  :

  Print method for `chain_structure_group`.

- [`print(`*`<cluster_choice>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.cluster_choice.md)
  : Print Method for cluster_choice

- [`print(`*`<mcml>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.mcml.md)
  : Print Method for mcml

- [`print(`*`<mcml_layer>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.mcml_layer.md)
  : Print Method for an mcml Layer

- [`print(`*`<mmm_compare>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.mmm_compare.md)
  : Print Method for mmm_compare

- [`print(`*`<nestimate_data>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.nestimate_data.md)
  : Print Method for nestimate_data

- [`print(`*`<net_association_rules>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_association_rules.md)
  : Print Method for net_association_rules

- [`print(`*`<net_bootstrap>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_bootstrap.md)
  : Print Method for net_bootstrap

- [`print(`*`<net_bootstrap_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_bootstrap_group.md)
  : Print Method for net_bootstrap_group

- [`print(`*`<net_cluster_diagnostics>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_cluster_diagnostics.md)
  : Print Method for net_cluster_diagnostics

- [`print(`*`<net_clustering>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_clustering.md)
  : Print Method for net_clustering

- [`print(`*`<net_gimme>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_gimme.md)
  : Print Method for net_gimme

- [`print(`*`<net_hon>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_hon.md)
  : Print Method for net_hon

- [`print(`*`<net_honem>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_honem.md)
  : Print Method for net_honem

- [`print(`*`<net_hypa>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_hypa.md)
  : Print Method for net_hypa

- [`print(`*`<net_link_prediction>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_link_prediction.md)
  : Print Method for net_link_prediction

- [`print(`*`<net_markov_order>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_markov_order.md)
  : Print Method for net_markov_order

- [`print(`*`<net_markov_order_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_markov_order_group.md)
  :

  Print method for `net_markov_order_group`

- [`print(`*`<net_markov_stability_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_markov_stability_group.md)
  :

  Print method for `net_markov_stability_group`

- [`print(`*`<net_mlvar>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_mlvar.md)
  : Print method for net_mlvar

- [`print(`*`<net_mmm>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_mmm.md)
  : Print Method for net_mmm

- [`print(`*`<net_mmm_clustering>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_mmm_clustering.md)
  : Print Method for MMM Clustering Attribute

- [`print(`*`<net_mogen>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_mogen.md)
  : Print Method for net_mogen

- [`print(`*`<net_mpt_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_mpt_group.md)
  :

  Print method for `net_mpt_group`

- [`print(`*`<net_nct>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_nct.md)
  : Print Method for net_nct

- [`print(`*`<net_path_dependence>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_path_dependence.md)
  :

  Print method for `net_path_dependence`

- [`print(`*`<net_permutation>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_permutation.md)
  : Print Method for net_permutation

- [`print(`*`<net_permutation_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_permutation_group.md)
  : Print Method for net_permutation_group

- [`print(`*`<net_reliability>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_reliability.md)
  : Print Method for net_reliability

- [`print(`*`<net_sequence_comparison>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_sequence_comparison.md)
  : Print Method for net_sequence_comparison

- [`print(`*`<net_stability>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_stability.md)
  : Print Method for net_stability

- [`print(`*`<net_stability_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_stability_group.md)
  : Print Method for net_stability_group

- [`print(`*`<net_transition_entropy>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_transition_entropy.md)
  :

  Print method for `net_transition_entropy`

- [`print(`*`<net_transition_entropy_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.net_transition_entropy_group.md)
  :

  Print method for `net_transition_entropy_group`

- [`print(`*`<netobject>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.netobject.md)
  : Print Method for Network Object

- [`print(`*`<netobject_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.netobject_group.md)
  : Print Method for Group Network Object

- [`print(`*`<netobject_ml>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.netobject_ml.md)
  : Print Method for Multilevel Network Object

- [`print(`*`<persistent_homology>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.persistent_homology.md)
  : Print persistent homology results

- [`print(`*`<q_analysis>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.q_analysis.md)
  : Print Q-analysis results

- [`print(`*`<simplicial_complex>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.simplicial_complex.md)
  : Print a simplicial complex

- [`print(`*`<summary.net_path_dependence>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.summary.net_path_dependence.md)
  :

  Print method for `summary.net_path_dependence`

- [`print(`*`<summary.net_transition_entropy>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.summary.net_transition_entropy.md)
  :

  Print method for `summary.net_transition_entropy`

- [`print(`*`<summary_chain_structure>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.summary_chain_structure.md)
  :

  Print method for `summary.chain_structure`.

- [`print(`*`<wtna_boot_mixed>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.wtna_boot_mixed.md)
  : Print Method for wtna_boot_mixed

- [`print(`*`<wtna_mixed>`*`)`](https://mohsaqr.github.io/Nestimate/reference/print.wtna_mixed.md)
  : Print Method for wtna_mixed

- [`print(`*`<state_freq>`*`)`](https://mohsaqr.github.io/Nestimate/reference/state_freq.md)
  [`plot(`*`<state_freq>`*`)`](https://mohsaqr.github.io/Nestimate/reference/state_freq.md)
  [`as.data.frame(`*`<state_freq>`*`)`](https://mohsaqr.github.io/Nestimate/reference/state_freq.md)
  : Print, Plot, and Convert a state_freq Object

- [`passage_time()`](https://mohsaqr.github.io/Nestimate/reference/passage_time.md)
  [`summary(`*`<net_mpt>`*`)`](https://mohsaqr.github.io/Nestimate/reference/passage_time.md)
  [`plot(`*`<net_mpt>`*`)`](https://mohsaqr.github.io/Nestimate/reference/passage_time.md)
  : Mean First Passage Times

- [`summary(`*`<boot_glasso>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.boot_glasso.md)
  : Summary Method for boot_glasso

- [`summary(`*`<chain_structure>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.chain_structure.md)
  :

  Tidy per-state summary of a `chain_structure`.

- [`summary(`*`<chain_structure_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.chain_structure_group.md)
  :

  Cross-group comparison of `chain_structure_group`.

- [`summary(`*`<cluster_choice>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.cluster_choice.md)
  : Summary Method for cluster_choice

- [`summary(`*`<mcml>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.mcml.md)
  : Summary Method for mcml

- [`summary(`*`<mmm_compare>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.mmm_compare.md)
  : Summary Method for mmm_compare

- [`summary(`*`<nest_initial_probs>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.nest_initial_probs.md)
  : Summary Method for Initial Probability Vectors

- [`summary(`*`<nest_transition_counts>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.nest_transition_counts.md)
  : Summary Method for Transition Count Matrices

- [`summary(`*`<nest_transition_matrix>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.nest_transition_matrix.md)
  : Summary Method for Transition Matrices

- [`summary(`*`<net_association_rules>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_association_rules.md)
  : Summary Method for net_association_rules

- [`summary(`*`<net_bootstrap>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_bootstrap.md)
  : Summary Method for net_bootstrap

- [`summary(`*`<net_bootstrap_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_bootstrap_group.md)
  : Summary Method for net_bootstrap_group

- [`summary(`*`<net_clustering>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_clustering.md)
  : Summary Method for net_clustering

- [`summary(`*`<net_gimme>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_gimme.md)
  : Summary Method for net_gimme

- [`summary(`*`<net_hon>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_hon.md)
  : Summary Method for net_hon

- [`summary(`*`<net_honem>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_honem.md)
  : Summary Method for net_honem

- [`summary(`*`<net_hypa>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_hypa.md)
  : Summary Method for net_hypa

- [`summary(`*`<net_link_prediction>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_link_prediction.md)
  : Summary Method for net_link_prediction

- [`summary(`*`<net_markov_order>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_markov_order.md)
  : Summary Method for net_markov_order

- [`summary(`*`<net_mlvar>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_mlvar.md)
  : Summary method for net_mlvar

- [`summary(`*`<net_mmm>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_mmm.md)
  : Summary Method for net_mmm

- [`summary(`*`<net_mogen>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_mogen.md)
  : Summary Method for net_mogen

- [`summary(`*`<net_nct>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_nct.md)
  : Summary Method for net_nct

- [`summary(`*`<net_path_dependence>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_path_dependence.md)
  :

  Summary method for `net_path_dependence`

- [`summary(`*`<net_permutation>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_permutation.md)
  : Summary Method for net_permutation

- [`summary(`*`<net_permutation_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_permutation_group.md)
  : Summary Method for net_permutation_group

- [`summary(`*`<net_reliability>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_reliability.md)
  : Summary Method for net_reliability

- [`summary(`*`<net_sequence_comparison>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_sequence_comparison.md)
  : Summary Method for net_sequence_comparison

- [`summary(`*`<net_stability>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_stability.md)
  : Summary Method for net_stability

- [`summary(`*`<net_stability_group>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_stability_group.md)
  : Summary Method for net_stability_group

- [`summary(`*`<net_transition_entropy>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.net_transition_entropy.md)
  :

  Summary method for `net_transition_entropy`

- [`summary(`*`<wtna_boot_mixed>`*`)`](https://mohsaqr.github.io/Nestimate/reference/summary.wtna_boot_mixed.md)
  : Summary Method for wtna_boot_mixed

- [`markov_stability()`](https://mohsaqr.github.io/Nestimate/reference/markov_stability.md)
  [`plot(`*`<net_markov_stability>`*`)`](https://mohsaqr.github.io/Nestimate/reference/markov_stability.md)
  : Markov Stability Analysis

- [`plot(`*`<boot_glasso>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.boot_glasso.md)
  : Plot Method for boot_glasso

- [`plot(`*`<chain_structure>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.chain_structure.md)
  :

  Plot method for `chain_structure`.

- [`plot(`*`<cluster_choice>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.cluster_choice.md)
  : Plot Method for cluster_choice

- [`plot(`*`<mmm_compare>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.mmm_compare.md)
  : Plot Method for mmm_compare

- [`plot(`*`<net_association_rules>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_association_rules.md)
  : Plot Method for net_association_rules

- [`plot(`*`<net_cluster_diagnostics>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_cluster_diagnostics.md)
  : Plot Method for net_cluster_diagnostics

- [`plot(`*`<net_clustering>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_clustering.md)
  : Plot Sequence Clustering Results

- [`plot(`*`<net_gimme>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_gimme.md)
  : Plot Method for net_gimme

- [`plot(`*`<net_honem>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_honem.md)
  : Plot Method for net_honem

- [`plot(`*`<net_markov_order>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_markov_order.md)
  : Plot Method for net_markov_order

- [`plot(`*`<net_mmm>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_mmm.md)
  : Plot Method for net_mmm

- [`plot(`*`<net_mmm_clustering>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_mmm_clustering.md)
  : Plot Method for MMM Clustering Attribute

- [`plot(`*`<net_mogen>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_mogen.md)
  : Plot Method for net_mogen

- [`plot(`*`<net_path_dependence>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_path_dependence.md)
  :

  Plot method for `net_path_dependence`

- [`plot(`*`<net_reliability>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_reliability.md)
  : Plot Method for net_reliability

- [`plot(`*`<net_sequence_comparison>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_sequence_comparison.md)
  : Plot Method for net_sequence_comparison

- [`plot(`*`<net_stability>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_stability.md)
  : Plot Method for net_stability

- [`plot(`*`<net_transition_entropy>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.net_transition_entropy.md)
  :

  Plot method for `net_transition_entropy`

- [`plot(`*`<persistent_homology>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.persistent_homology.md)
  : Plot Persistent Homology

- [`plot(`*`<q_analysis>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.q_analysis.md)
  : Plot Q-Analysis

- [`plot(`*`<simplicial_complex>`*`)`](https://mohsaqr.github.io/Nestimate/reference/plot.simplicial_complex.md)
  : Plot a Simplicial Complex

- [`plot_mosaic()`](https://mohsaqr.github.io/Nestimate/reference/plot_mosaic.md)
  : Draw a Marimekko / Mosaic Plot from a Tidy Data Frame

- [`plot_state_frequencies()`](https://mohsaqr.github.io/Nestimate/reference/plot_state_frequencies.md)
  : Plot State Frequency Distributions
