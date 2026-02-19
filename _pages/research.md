---
layout: archive
title: "Research"
permalink: /research/
author_profile: true
redirect_from:
  - /research
---


## Working Papers

***"Cross-Fitting-Free Debiased Machine Learning with Multiway Dependence"***  \
with **Harold D. Chiang** \
[arXiv](https://arxiv.org/abs/2602.11333)
<details>
<summary>Abstract</summary>
This paper develops an asymptotic theory for two-step debiased machine learning (DML) estimators in generalised method of moments (GMM) models with general multiway clustered dependence, without relying on cross-fitting. While cross-fitting is commonly employed, it can be statistically inefficient and computationally burdensome when first-stage learners are complex and the effective sample size is governed by the number of independent clusters. We show that valid inference can be achieved without sample splitting by combining Neyman-orthogonal moment conditions with a localisation-based empirical process approach, allowing for an arbitrary number of clustering dimensions. The resulting DML-GMM estimators are shown to be asymptotically linear and asymptotically normal under multiway clustered dependence. A central technical contribution of the paper is the derivation of novel global and local maximal inequalities for general classes of functions of sums of separately exchangeable arrays, which underpin our theoretical arguments and are of independent interest.
</details>


***"Inference in High-Dimensional Panel Models: Two-Way Dependence and Unobserved Heterogeneity"***  \
Best Student Paper Award at *2024 Midwest Econometrics Group Meeting* \
[arXiv](https://arxiv.org/abs/2504.18772) &nbsp; [latest version](https://kaichengchen.github.io/TW_DML_LASSO_CRE.pdf) &nbsp; [code & replication](http://kaichengchen.github.io/twlasso_paneldml_replication.zip)
<details>
<summary>Abstract</summary>
Panel data allows for the modeling of unobserved heterogeneity, significantly raising the number of nuisance parameters and making high dimensionality a practical issue. Meanwhile, temporal and cross-sectional dependence in panel data further complicates high-dimensional estimation and inference. This paper proposes a toolkit for high-dimensional panel models with large cross-sectional and time sample sizes. To reduce the dimensionality, I propose a variant of LASSO for two-way clustered panels. While being consistent, the convergence rate of LASSO is slow due to the cluster dependence, rendering inference challenging in general. Nevertheless, asymptotic normality can be established in a semiparametric moment-restriction model by leveraging a clustered-panel cross-fitting approach and, as a special case, in a partial linear model using the full sample. In an exercise of estimating multiplier using panel data, I demonstrate how high dimensionality could be hidden and the proposed toolkit enables flexible modeling and robust inference.
</details>

***"Identification of Partial Effects with Endogenous Controls"***\
with **Kyoo il Kim** \
[arXiv](https://arxiv.org/abs/2401.14395) &nbsp; [latest version](https://kaichengchen.github.io/endogenous_control.pdf)
 <details>
<summary>Abstract</summary>
Control variables are routinely treated as exogenous, yet in many empirical settings they are themselves endogenous. This creates a dilemma: omitting controls can leave the treatment endogenous, while including them may contaminate identification. The problem is not resolved with an instrumental variable when it is only conditionally valid. We show that average responses to the treatment remain identified under a rank condition—measurable separability—that accommodates endogenous controls. For parametric models, our approach amounts to estimating a nonparametric model that nests the parametric specification. For nonparametric models, our results imply that endogenous controls are generally innocuous under standard identification conditions, except in the presence of ``bad controls''. We further propose a test for endogenous controls. Simulation results and an empirical application demonstrate this prevalent issue and provide practical implications of our methods.
 </details>


***"Another Look at the Linear Probability Model and Nonlinear Index Models"***\
 with **Robert S. Martin** and **Jeffrey M. Wooldridge** (Revise & Resubmit at *Econometric Reviews*) \
[arXiv](https://arxiv.org/abs/2308.15338) &nbsp; [latest version](https://kaichengchen.github.io/LPM_CMW.pdf) &nbsp; [code & replication](https://kaichengchen.github.io/lpm_simulation_post.rar)
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


