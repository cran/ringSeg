# ringSeg 0.1.0

* Initial release.
* RING rank-based, asymptotic distribution-free change-point detection
  (Zhou & Chen, 2025, IEEE Transactions on Information Theory,
  <doi:10.1109/TIT.2025.3575858>): scans a sequence for a single change-point
  or a changed interval using the weighted (WR), max-type (MR), and generalized
  (TR) rank statistics, with analytic distribution-free p-values, an optional
  skewness correction, and optional permutation p-values.
* Main functions: `ring_cpd()` (data / distance in), `rcpd()` (rank matrix in),
  `ring_graph()` (build the RING k-NN rank graph), and `Rise_Rank()`.
