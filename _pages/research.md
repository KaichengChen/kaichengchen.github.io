---
layout: archive
title: "Research"
permalink: /research/
author_profile: true
redirect_from:
  - /research
---


## Working Papers

***"Inference in High-Dimensional Panel Models: Two-Way Dependence and Unobserved Heterogeneity"***  \
**Best Student Paper Award**, 2024 Midwest Econometrics Group Meeting \
[arXiv](https://arxiv.org/abs/2504.18772) &nbsp; [latest version](https://kaichengchen.github.io/TW_DML_LASSO_CRE.pdf) 
<details>
<summary>Abstract</summary>
Panel data allows for the modeling of unobserved heterogeneity, which significantly increases the number of nuisance parameters, making high dimensionality a practical issue rather than just a theoretical concern. However, unobserved heterogeneity, along with potential temporal and cross-sectional dependence in panel data, further complicates estimation and inference for high-dimensional models. This paper proposes a toolkit for robust estimation and inference in high-dimensional panel models with large cross-sectional and time sample sizes. To reduce the dimensionality, I propose a weighted LASSO using two-way cluster-robust penalty weights. Due to the cluster dependence, the rate of convergence is slow even in an oracle case. Nevertheless, by leveraging a clustered-panel cross-fitting approach for bias correction, asymptotic normality can be established for the low-dimensional vector of the estimated parameters. As a special case, inferential theories are also established using the full sample in a partial linear model with unobserved time and unit effects. In a panel estimation of the government spending multiplier, I demonstrate how high dimensionality can be hidden and how the proposed toolkit enables flexible modeling and robust inference.
</details>

***"Identification of Partial Effects with Endogenous Controls"***\
with **Kyoo il Kim** \
[arXiv](https://arxiv.org/abs/2401.14395) &nbsp; [latest version](https://kaichengchen.github.io/endogenous_control.pdf)
 <details>
<summary>Abstract</summary>
Exogeneity of the treatment needed for identification are often achieved by conditioning. While control variables are explicitly or implicitly assumed to be exogenous, it is common to encounter endogenous controls in practice. It brings a dilemma: without controlling, the treatment may be endogenous; with controlling, the endogeneity of controls may pollute the identification. The problem is not solved with an instrumental variable when it is only conditionally valid and controls are endogenous. We provide identification results for local average response under an extra measurable separability condition between the treatment and the controls. Noticeably, this condition permits the controls to be dependent on the treatment. The results apply to a wide class of models ranging from linear to non-separable ones. Monte Carlo simulations exemplify this prevalent issue and demonstrate the performance of the proposed methods in finite sample.
 </details>


***"Another Look at the Linear Probability Model and Nonlinear Index Models"***\
 with **Robert S. Martin** and **Jeffrey M. Wooldridge**\
 [arXiv](https://arxiv.org/abs/2308.15338) &nbsp; [latest version](https://kaichengchen.github.io/LPM_CMW.pdf)
 <details>
<summary>Abstract</summary>
We reassess the use of linear models for binary responses, focusing on average partial effects (APEs). We confirm that under certain conditions, linear projection parameters correspond to APEs even when the true model is nonlinear. Simulations demonstrate a large fraction of fitted values in [0, 1] is neither necessary nor sufficient for OLS to approximate the APEs. To reduce bias, excluding observations with fitted values outside [0, 1] has been proposed. We show that iteratively trimming the sample is equivalent to nonlinear least squares estimation of a piece-wise linear (ramp) model, for which we establish consistency and asymptotic normality results.
</details>


## Publication
***"Fixed-b Asymptotics for Panel Models with Two-Way Clustering"***\
with **Timothy J. Vogelsang** (*Journal of Econometrics*, 2024, 244(1): 105831) \
[arXiv](https://arxiv.org/abs/2309.08707) &nbsp; [published version](https://www.sciencedirect.com/science/article/abs/pii/S0304407624001763)\
STATA command: ``xtregtfb``. Installation: type ``net from https://kaichengchen.github.io/statafile/`` in STATA 
<details>
<summary>Abstract</summary>
This paper studies a cluster robust variance estimator proposed by Chiang, Hansen and Sasaki (2024) for linear panels. First, we show algebraically that this variance estimator (CHS estimator, hereafter) is a linear combination of three common variance estimators: the one-way unit cluster estimator, the "HAC of averages" estimator, and the
"average of HACs" estimator. Based on this finding, we obtain a fixed-b asymptotic result for the CHS estimator and corresponding test statistics as the cross-section and time sample sizes jointly go to infinity. Furthermore, we propose two simple bias-corrected versions of the variance estimator and derive the fixed-b limits. In a
simulation study, we find that the two bias-corrected variance estimators along with fixed-b critical values provide improvements in finite sample coverage probabilities. We illustrate the impact of bias-correction and use of the fixed-b critical values on inference in an empirical example on the relationship between industry profitability and market concentration.
</details>


