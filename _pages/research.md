---
layout: archive
title: "Research"
permalink: /research/
author_profile: true
redirect_from:
  - /research
---

## Job Market Paper

***"Inference in High-Dimensional Panel Models: Two-Way Dependence and Unobserved Heterogeneity"***  
**Best Student Paper Award**, 2024 Midwest Econometrics Group Meeting \
[latest version](https://kaichengchen.github.io/TW_DML_LASSO_CRE.pdf)
<details>
<summary>Abstract</summary>
Panel data allows for the modeling of unobserved heterogeneity, which significantly increases the number of nuisance parameters, making high dimensionality a practical issue rather than just a theoretical concern. However, unobserved heterogeneity, along with potential temporal and cross-sectional dependence in panel data, further complicates estimation and inference for high-dimensional models. This paper proposes a toolkit for robust estimation and inference in high-dimensional panel models with large cross-sectional and time sample sizes. To reduce the dimensionality, I propose a weighted LASSO using two-way cluster-robust penalty weights. Due to the cluster dependence driven by the underlying components, the rate of convergence is slow even in an oracle case. Nevertheless, by leveraging a clustered-panel cross-fitting approach for bias-correction, the asymptotic normality on low-dimensional parameters can be established using the weighted LASSO for nuisance estimation. As a special case, in a partial linear model with non-additive unobserved time and unit effects, inferential results are also established using the full sample. In a panel estimation of the government spending multiplier, I demonstrate how high dimensionality can be hidden and how the proposed toolkit enables flexible modeling and robust inference.
</details>


## Working Papers

**["Identification of Nonseparable Models with Endogenous Control Variables"](https://arxiv.org/abs/2401.14395)**\
with **Kyoo il Kim**
 <details>
<summary>Abstract</summary>
Identification of partial effects relies on some exogeneity conditions of the targeted treatment, which is often achieved through including relevant control variables. While these controls are implicitly or explicitly assumed to be exogenous, it is common to encounter endogenous control variables in practice. It brings a dilemma: without controlling, both the unobserved determinants of the outcome and the relevant controls cause the endogeneity issue for the treatment; with controlling, the endogeneity of controls will pollute the identification even with the conditional independence. Either way, due to the lack of identification, estimations assuming exogeneity of the controls are biased and the inference are rendered invalid. The problem is not solved with an instrumental variable when the IV is only conditionally valid and controls are endogenous. We provide an alternative identification for both cases under an extra measurable separability condition between the treatment and the controls. Noticeably, this condition permits the controls to be influenced by the treatment, effectively allowing for some types of bad controls. The results apply to a wide class of models including linear, nonlinear, and non-separable models. Monte Carlo simulations exemplify the bias of estimations based on exogeneity assumption on the control when they are actually endogenous, and the proposed identification methods paired with usual nonparametric estimators perform well in finite sample. We revisit empirical studies published in top economics journals and show how our methods matter in practice.
</details>


**["Another Look at the Linear Probability Model and Nonlinear Index Models"](https://arxiv.org/abs/2308.15338)**\
 with **Robert S. Martin** and **Jeffrey M. Wooldridge**
 <details>
<summary>Abstract</summary>
We reassess the use of linear models for binary responses, focusing on average partial effects (APEs). We confirm that under certain conditions, linear projection parameters correspond to APEs even when the true model is nonlinear. Simulations demonstrate a large fraction of fitted values in [0, 1] is neither necessary nor sufficient for OLS to approximate the APEs. To reduce bias, excluding observations with fitted values outside [0, 1] has been proposed. We show that iteratively trimming the sample is equivalent to nonlinear least squares estimation of a piece-wise linear (ramp) model, for which we establish consistency and asymptotic normality results.
</details>


## Publication
**["Fixed-b Asymptotics for Panel Models with Two-Way Clustering"](https://arxiv.org/abs/2309.08707)**\
with **Timothy J. Vogelsang** (*Journal of Econometrics*, 2024, 244(1): 105831) \
STATA command: ``xtregtfb``. Installation: type ``net from https://kaichengchen.github.io/statafile/`` in STATA 
<details>
<summary>Abstract</summary>
This paper studies a cluster robust variance estimator proposed by Chiang, Hansen and Sasaki (2024) for linear panels. First, we show algebraically that this variance estimator (CHS estimator, hereafter) is a linear combination of three common variance estimators: the one-way unit cluster estimator, the "HAC of averages" estimator, and the
"average of HACs" estimator. Based on this finding, we obtain a fixed-b asymptotic result for the CHS estimator and corresponding test statistics as the cross-section and time sample sizes jointly go to infinity. Furthermore, we propose two simple bias-corrected versions of the variance estimator and derive the fixed-b limits. In a
simulation study, we find that the two bias-corrected variance estimators along with fixed-b critical values provide improvements in finite sample coverage probabilities. We illustrate the impact of bias-correction and use of the fixed-b critical values on inference in an empirical example on the relationship between industry profitability and market concentration.
</details>


