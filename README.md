# Poisson3armNI: R codes for "New Approaches for Testing Non-inferiority for Three-arm Trials with Poisson Distributed Outcomes"

A legacy version of the software code exists in the following webpage: 

https://github.com/Shrabanti87/PoisNI3arm   

The codes are provided to generate sample size and power for both Frequentist and Bayesian approaches. The following notations are used in this package:

•	lamE (or l_E), l_R (or lamR), and l_P (or lamP): The Poisson rate parameters for the arm E, R, and P, respectively

•	nP: The sample size in the placebo arm (P)

•	nR: The sample size in the reference arm (R)

•	nE: The sample size in the experimental arm (E)

•	theta: effect retention parameter

•	r.alloc (or alloc): Allocation vector which can be (1:1:1), (2:2:1) or (3:2:1) for nE:nR:nP

•	n: Total sample size

We give brief description of the R files below:

1. freq_power_sample

This function calculates the Frequentist sample size for a given value of theta, allocation, alpha, l_P, l_R, and l_E.

Arguments

n_total: Maximum number of n to get the sample size

alloc: Allocation vector

l_E: lambda of the arm E

l_R: lambda of the arm R

l_P: lambda of the arm P

theta: effect retention parameter

Output: Power, power curves for different allocations, power curves for different theta, and minimum sample size of the arm P satisfying power>=1-beta 


2. type1error_fullybayesian.R

This function calculates the estimated type I error under exact Bayesian approach for different values of theta, calculated sample size, allocation, and l_E.

Arguments

n: sample size

alloc: Allocation vector

l_E: lambda of the arm E

l_R: lambda of the arm R

l_P: lambda of the arm P

theta: effect retention parameter

aE, bE, aR, bR, aP, bP: Bayesian gamma prior parameters for three arms

Output: Estimated type I error for exact Bayesian approach

3. power_poisson_functions.R

This file contains the functions that need to be sourced first to calculate sample size under Frequentist and Bayesian approaches. The main functions, arguments and outputs are described below:

a) power_freq

This function calculates the Frquentist power for a given sample size, allocation, lamE and theta.

Arguments

n: Sample size in the arm P

r.alloc: Allocation vector

lamE: lambda of the arm E

theta: effect retention parameter

Output

power: Power for the specified values of the parameters in the arguments

b) samplesize_fn_freq

This function calculates the Frequentist sample size for a particular value of theta for a given allocation and lamE.

Arguments

r.alloc: Allocation vector

lamE: lambda of the arm E

theta: effect retention parameter

Output

n: Minimum sample size of the arm P satisfying power>=0.8

c) power_fbayes

This function calculates the fully Bayesian power for a given sample size, allocation, lamE and theta.

Arguments

n: Sample size in the arm P

r.alloc: Allocation vector

lamE: lambda of the arm E

theta: effect retention parameter

Output

power: Power for the specified values of the parameters in the arguments

d) samplesize_fn_fbayes

This function calculates the fully Bayesian sample size for a particular value of theta for a given allocation and lamE.

Arguments

r.alloc: Allocation vector

lamE: lambda of the arm E

theta: effect retention parameter

a_max: maximum range of search for the Bayesian sample size

Output

n: Minimum sample size of the arm P satisfying power>=0.8

e) power_approxbayes

This function calculates the approximation-based Bayesian power for a given sample size, allocation, lamE and theta.

Arguments

n: Sample size in the arm P

r.alloc: Allocation vector

lamE: lambda of the arm E

theta: effect retention parameter

Output

power: Power for the specified values of the parameters in the arguments

f) samplesize_fn_approxbayes

This function calculates the sample size under approximation-based Bayesian approach for a particular value of theta for a given allocation and lamE.

Arguments

r.alloc: Allocation vector

lamE: lambda of the arm E

theta: effect retention parameter

Output

n: Minimum sample size of the arm P satisfying power>=0.8

4. samplesize_calc.R
This file calculates the sample size for 80% power under the Frequentist, fully Bayesian and approximation-based Bayesian approaches across different values of lamE, theta and allocations.

5. freq.simulatedpower.R
This file calculates the Frequentist simulated power across different values of lamE, theta and allocations using the calculated sample size.

6. full.bayes.simulatedpower.R
This file calculates the fully Bayesian simulated power across different values of lamE, theta and allocations using the calculated sample size.

7. approx.bayes.simulatedpower.R
This file calculates the simulated power under approximation-based Bayesian approach across different values of lamE, theta and allocations using the calculated sample size.
