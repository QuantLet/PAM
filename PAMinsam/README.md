<div style="margin: 0; padding: 0; text-align: center; border: none;">
<a href="https://quantlet.com" target="_blank" style="text-decoration: none; border: none;">
<img src="https://github.com/StefanGam/test-repo/blob/main/quantlet_design.png?raw=true" alt="Header Image" width="100%" style="margin: 0; padding: 0; display: block; border: none;" />
</a>
</div>

```
Name of Quantlet: PAMinsam

Published in: submitted to N/A

Description: ‘Performs the Penalized Adaptive Method (PAM), a combination of propagation-separation approach and a SCAD penalty, to fit a model with possible structural changes on a given dataset. Compares the fit to models of Cochrane and Piazzesi (2005) and Ludvigson and Ng (2009). The input data are monthly observations of k-year excess bond risk premia, k = 2, 3, 4, 5, forward rates and other pre-defined macro variables. Computes RMSE, MAE, R^2 and adjusted R^2 for the fitted models. Plots the time series of the fitted vs. observed values of bond risk excess premia.'

Keywords: ‘linear model, regression, SCAD penalty, bic, time varying, change point, bootstrap, plot, visualization, nonstationary, financial, returns’

See also: ‘PAMsimLR, PAMsimCP, PAMCocPia, PAMoutsam’

Author: Lenka Zboňáková

Submitted: 23 May 2018 by Lenka Zboňáková

Datafile: BRP_data.csv, LN_macrodata_transformed.csv

Input: 

- start.date: Starting year of observations for PAM and CP

- end.date: Ending year of observations for PAM and CP

- start.dateLN: Starting year of observations for LN

- end.dateLN: Ending year of observations for LN

- Y.real: Response variable (k-year excess bond risk premia, k = 2, 3, 4, 5)

- n.years: Number of years as increment between successive subintervals

- n.boot: Number of bootstrap loops

- mb.type: Distribution of multipliers (Bound, Exp or Pois)

```
<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/PAM/master/PAMinsam/PAMinsam-1.png" alt="Image" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/PAM/master/PAMinsam/PAMinsam-2.png" alt="Image" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/PAM/master/PAMinsam/PAMinsam-3.png" alt="Image" />
</div>

