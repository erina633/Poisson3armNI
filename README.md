# Poisson3armNI
R codes for " New Approaches for Testing Non-inferiority for Three-arm Trials with Poisson Distributed Outcomes"
PoisNI3arm
R codes for " New Approaches for Testing Non-inferiority for Three-arm Trials with Poisson Distributed Outcomes"
The codes are provided to generate sample size and power for both Frequentist and Bayesian approaches. The following notations are used in this package:
•	lamE (or l_E), lamR (or l_R),  and lamP (or l_P): The Poisson rate parameter in the arm E, R, and P respectively
•	nP: The sample size in the placebo arm (P)
•	nR: The sample size in the reference arm (R)
•	nE: The sample size in the experimental arm (E)
•	theta: effect retention parameter
•	r.alloc (or alloc): Allocation vector which can be (1:1:1), (2:2:1) or (3:2:1) for nE:nR:nP
•	N: Total sample size
We give brief description of the R files below:

1. type1error_fullybayesian.R
This function calculates the estimated type I error under exact Bayesian approach for different values of theta, calculated sample size, allocation, and lamE.
Arguments
•	n: sample size
•	r.alloc: Allocation vector
•	lamE: lambda of the arm E
•	theta: effect retention parameter
Output
Estimated type I error
