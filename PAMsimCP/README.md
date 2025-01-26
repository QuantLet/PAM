<div style="margin: 0; padding: 0; text-align: center; border: none;">
<a href="https://quantlet.com" target="_blank" style="text-decoration: none; border: none;">
<img src="https://github.com/StefanGam/test-repo/blob/main/quantlet_design.png?raw=true" alt="Header Image" width="100%" style="margin: 0; padding: 0; display: block; border: none;" />
</a>
</div>

```
Name of Quantlet: PAMsimCP

Published in: submitted to N/A

Description: ‘Performs the Penalized Adaptive Method (PAM), a combination of propagation-separation approach and a SCAD penalty, to fit a model to a simulated data with a pre-defined change point. Computes the percentage of correctly identified change points over a specified number of scenarios.’

Keywords: ‘linear model, regression, SCAD penalty, bic, change point, simulations, bootstrap’

See also: ‘PAMsimLR, PAMCocPia, PAMinsam, PAMoutsam’

Author: Lenka Zboňáková

Submitted: 23 May 2018 by Lenka Zboňáková

Input: 
- n.obs   : Number of observations
- n.par   : Number of parameters
- n.sim   : Number of simulated scenarios
- r       : Correlation parameter of the design matrix X
- sd.eps  : Standard deviation of the error term
- cp1.seq : Sequence of change points
- n.boot  : Number of bootstrap loops
- K       : Sequence of increments between adjacent subintervals
- mb.type : Distribution of multipliers (Bound, Exp or Pois)

```
