# Exploring deterministic pointsets in computational Bayesian methodology

Applying Bayesian methods to models with complex hierarchical structures can be computationally expensive. Integrated nested Laplace approximations (INLA) is a Bayesian method that is known to be more efficient than other Bayesian methods. But INLA employs deterministic point sets such as grid structures, empirical bayes (EB) and central composite design (CCD) to compute posterior estimates. This project explores alternative deterministic point sets within INLA, specifically point sets with low discrepancy, and evaluate their perfomance on a range of models with different complexities. 

In comparison to INLA's default exploration schemes (grid, EB, CCD), do the properties of Low Discrepancy Sequences (LDS) improve performance when applied to Bayesian models with complex hierarchical structures?
