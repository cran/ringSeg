# ringSeg 0.1.1

* Preserve very small uncorrected max-type p-values when combining tails.
* Return missing analytic p-values below the large-threshold tail range.
* Stabilize scans with nearly constant edge weights and exclude invalid
  saddlepoints from skewness-corrected integrals.
* Stabilize the saddlepoint calculation at zero and near-zero skewness.
* Use the uncorrected analytic component when skewness extrapolation has no
  usable anchors or fails numerically.
* Share the rank weights of tied neighbours that straddle the k-NN cutoff.
* Evaluate both components of the corrected max-type p-value at the max-type
  scan threshold, for both scan alternatives.
* Restrict changed-interval maxima to eligible intervals. Handle zero null
  variance consistently in observed and permuted scans. Return missing
  analytic p-values for fixed bounds or zero null variance.
* Rank distances directly to avoid creating ties through floating-point
  subtraction.

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
