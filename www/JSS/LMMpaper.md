LMM paper - Bates, Maechler, Bolker and Walker

1. Introduction (done)
  * Motivation
  * Model and notation

2. Derivation of the modified Henderson's MME and profiled log-likelihood (done)

3. Implementation of PLS solution in pure R (and in Julia)
  * lme4-1.0 uses compiled code
  * getME function to extract objects of interest
	* VarCorr, pre-computation of Z'WZ, ...

4. Special cases
  * a single scalar r.e. term -> diagonal Lambda and L
	* a single vector-valued r.e. term -> block-diagonal Lambda and L
	* nested grouping factors -> no fill-in, permutation on post-order only
	* scalar, not necessarily nested model -> diagonal updates

5. Profiling parameters and graphical presentation

6. Modularization?

