{smcl}
{* *! version 1.0 Apr 17 2024}{...}
{title:Title}

{phang}
{bf:xtregtfb} {hline 2} provides bias-correction through fixed-b approximation for the two-way clustering robust standard error.

{marker syntax}{...}
{title:Syntax}
 
{p 4 17 2}
{cmd:xtregtfb}
{it:depvar}
[{it:indepvars}]
{ifin}
[{cmd:,} 
{bf:{ul:noc}onstant} {bf:fe} 
{cmd:se(}{it:setype}{cmd:)}
{cmd:lag(}{it:bandwidth}{cmd:)}
{cmd:level(}{it:confidencetwoside}{cmd:)}
{cmd:bm(}{it:bmincrement}{cmd:)}
{cmd:rep(}{it:repforsimcv}{cmd:)}
{cmd:whichvar(}{it:#}{cmd:)}]


{marker description}{...}
{title:Description}

{phang}
{cmd:xtregtfb} executes various bias-correction approaches to the CHS variance estimator based on {browse "https://arxiv.org/abs/2309.08707": Chen and Vogelsang (2023)}. In default, the bias-corrected CHS (BCCHS) standard error is reported. The simulated fixed-b critical values are reported only if {bf:se(string)} is specified with {bf:bcchsfb} or {bf:dkafb}, which can slow down the computation. More details and other options are described below.



{marker options}{...}
{title:Options}

{phang}
{bf:{ul:noc}onstant} {space 8}suppress constant term.
{p_end}

{phang}
{bf:fe} {space 16}fixed effects by within-transformation.
{p_end}

{phang}
{bf:se(string)} {space 8}specify the variance estimator reported in {bf:e(V)} and used for constructing the reported t-statistics: can be {bf:chs}, {bf:bcchs} (default), {bf:dka}, {bf:bcchsfb}, {bf:dkafb}.
{p_end}

{phang}
{bf:lag(integer)} {space 6}choose the bandwidth (any integer between 1 and the time sample size; by definition, 1 corresponds to CGM s.e.) for the Bartlett kernel; the default bandwidth is reported as {bf:e(Mhat)} is chosen using Eq(6.2) and Eq(6.4) of Andrews(1991) assuming AR(1) process, with 0 weight given to the constant term and other weights equal to the inverse squared variances of the estimated AR(1) processes. 
{marker options}{...}
{p_end}

{phang}
{bf:level(real)} {space 7}choose the significance level (of a two-sided test) for which the corresponding simulated fixed-b critical values are reported, if {bf:bcchsfb} or {bf:dkafb} is specified in {bf:se(string)}; the default is 0.05.
{p_end}

{phang}
{bf:whichvar(integer)} {space 1}specify the variable for which the simulated fixed-b critical values will be reported: the critical values are regressor-specific and the option is the sequence order among indepvar; the default is 1 (the first non-constant regressor).
{p_end}

{phang}
{bf:rep(integer)} {space 4}specify the number of replications for simulating the fixed-b critical values; the default is 2000.
{p_end}

{phang}
{bf:bm(integer)} {space 7}specify the number of increment for approximating the standard Wiener process; the default is 1000.
{p_end}


{title:Reference}

{p 4 8}Fixed-b Asymptotics for Panel Models with Two-Way Clustering. {browse "https://arxiv.org/abs/2309.08707": Chen and Vogelsang (2023)}.
{p_end}

