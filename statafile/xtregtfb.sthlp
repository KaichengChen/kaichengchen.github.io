{smcl}
{* *! version 1.0 Apr 12 2024}{...}
{title:Title}

{phang}
{bf:xtregtfb} {hline 2} provides bias-correction through fixed-b approximation for the two-way clustering robust standard error proposed by {browse "https://arxiv.org/abs/2201.11304":Chiang, Hansen, and Sasaki (2022)}.


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
{cmd:cvsim(}{it:simrep}{cmd:)}
{cmd:whichvar(}{it:#}{cmd:)}]


{marker description}{...}
{title:Description}

{phang}
{cmd:xtregtfb} executes bias-correction approaches to the CHS variance estimator based on {browse "https://arxiv.org/abs/2309.08707": Chen and Vogelsang (2023)}.



{marker options}{...}
{title:Options}

{phang}
{bf:{ul:noc}onstant} {space 8}suppress constant term.
{p_end}

{phang}
{bf:fe} {space 16}fixed effects by within-transformation.
{p_end}

{phang}
{bf:se(string)} {space 8}specify the variance estimator reported in e(V) and used for constructing the reported t-statistics: can be chs, bchs (default), dka, bchsfb, dkafb.
{p_end}

{phang}
{bf:lag(integer)} {space 6}choose the bandwidth (any integer between 1 and the time sample size; by definition, 1 corresponds to CGM s.e.) for the Bartlett kernel; the default bandwidth is chosen using Eq(6.2) and Eq(6.4) of Andrews(1991) assuming AR(1) process, with 0 weight given to the constant term and other weights equal to the inverse squared variances of the estimated AR(1) processes. 
{marker options}{...}
{p_end}

{phang}
{bf:level(real)} {space 7}choose the significance level (of a two-sided test) for which the corresponding simulated fixed-b critical values are reported, if bchsfb or dkafb is specified in se(); the default is 0.05.
{p_end}

{phang}
{bf:whichvar(integer)} {space 1}specify the variable for which the simulated fixed-b critical values will be reported: the critical values are regressor-specific and the option is the sequence order among indepvar; the default is 1 (the first non-constant regressor).
{p_end}

{phang}
{bf:cvsim(integer)} {space 4}specify the number of replications for simulating the fixed-b critical values; the default is 2000.
{p_end}

{phang}
{bf:bm(integer)} {space 7}specify the number of increment for approximating the standard Wiener process; the default is 1000.
{p_end}


{title:Reference}

{p 4 8}Fixed-b Asymptotics for Panel Models with Two-Way Clustering. {browse "https://arxiv.org/abs/2309.08707": Chen and Vogelsang (2023)}.
{p_end}

