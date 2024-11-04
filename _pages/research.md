---
layout: archive
title: "Research"
permalink: /research/
author_profile: true
redirect_from:
  - /research
---

## Job Market Paper

**["Inference in High-Dimensional Panel Models: Two-Way Dependence and Unobserved Heterogeneity"](https://kaichengchen.github.io/TW_DML_LASSO_CRE.pdf)**  
**Best Student Paper Award**, 2024 Midwest Econometrics Group Meeting
<details>
<summary>Abstract</summary>
Panel data allows for the modeling of unobserved heterogeneity, which significantly increases the number of nuisance parameters, making high dimensionality a practical issue rather than just a theoretical concern. However, unobserved heterogeneity, along with potential two-way dependence in panel data, further complicates estimation and inference for high-dimensional models. This paper proposes a toolkit for robust estimation and inference in high-dimensional panel models with large cross-sectional and time sample sizes. For estimation, I propose a weighted LASSO method with a theoretically driven penalty level and weights. Due to the high dimensionality and the two-way cluster dependence driven by the underlying components, the rate of convergence is slow even in an oracle case. Nonetheless, by leveraging a cross-fitting approach that is robust to panel data dependence, I show that it is possible to establish inferential theory on low-dimensional treatment parameters using the weighted LASSO for nuisance estimation. Additionally, I address the challenges posed by unobserved heterogeneity, which introduces a subtle issue for cross-fitting. Strategies and implications on the sparsity condition under various scenarios are discussed. In a panel estimation of the government spending multiplier, I demonstrate how high dimensionality can be hidden and how the proposed toolkit enables flexible modeling and robust inference.
</details>


## Working Papers

**["Identification of Nonseparable Models with Endogenous Control Variables"](https://arxiv.org/abs/2401.14395)**\
with **Kyoo il Kim**
 <details>
<summary>Abstract</summary>
We study identification of the treatment effects in a class of nonseparable models with the presence of potentially endogenous control variables. We show that given the treatment variable and the controls are measurably separated, the usual conditional independence condition or availability of excluded instrument suffices for identification.
</details>


**["Another Look at the Linear Probability Model and Nonlinear Index Models"](https://arxiv.org/abs/2308.15338)**\
 with **Robert S. Martin** and **Jeffrey M. Wooldridge**
 <details>
<summary>Abstract</summary>
We reassess the use of linear models to approximate response probabilities of binary outcomes, focusing on average partial effects (APEs). We confirm that linear projection parameters coincide with APEs in certain scenarios. Through simulations, we identify other cases where the linear projection does or does not approximate APEs and find that having a large fraction of fitted values in [0,1] is neither necessary nor sufficient. We also show nonlinear least squares estimation of the ramp model is consistent and asymptotically normal and is equivalent to using OLS on an iteratively trimmed sample to reduce bias. Our findings offer practical guidance for empirical research.
</details>


## Publication
**["Fixed-b Asymptotics for Panel Models with Two-Way Clustering"](https://urldefense.com/v3/__https://kwnsfk27.r.eu-west-1.awstrack.me/L0/https:*2F*2Fauthors.elsevier.com*2Fc*2F1jeqY15DjiIwZZ/1/01020191818f3bc6-a9f81387-b4b3-482a-a828-0707e93ed2c4-000000/wgqtxAAfCiwpXf66-aTMLrirVwk=388__;JSUlJQ!!HXCxUKc!xTVZ9jtRRyLUoZNF9HpSpbPWnYUM1OJAeOIIOobZushz2B02iIcvTE4gUKEqz_JnaVF0mJoy0PN24OqVcREH$)**\
with **Timothy J. Vogelsang** (*Journal of Econometrics*, 2024, 244(1): 105831) \
STATA command: ``xtregtfb``. Installation: type ``net from https://kaichengchen.github.io/statafile/`` in STATA 
<details>
<summary>Abstract</summary>
This paper studies a cluster robust variance estimator proposed by Chiang, Hansen and Sasaki (2024) for linear panels. First, we show algebraically that this variance estimator (CHS estimator, hereafter) is a linear combination of three common variance estimators: the one-way unit cluster estimator, the "HAC of averages" estimator, and the
"average of HACs" estimator. Based on this finding, we obtain a fixed-b asymptotic result for the CHS estimator and corresponding test statistics as the cross-section and time sample sizes jointly go to infinity. Furthermore, we propose two simple bias-corrected versions of the variance estimator and derive the fixed-b limits. In a
simulation study, we find that the two bias-corrected variance estimators along with fixed-b critical values provide improvements in finite sample coverage probabilities. We illustrate the impact of bias-correction and use of the fixed-b critical values on inference in an empirical example on the relationship between industry profitability and market concentration.
</details>


