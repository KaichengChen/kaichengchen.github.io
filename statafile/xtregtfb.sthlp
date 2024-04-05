{smcl}
{* *! version 1.1.1}{...}
{title:Title}

{phang}
{bf:xtregtfb} {hline 2} provides fixed-b adjustment/bias-correction for the two-way clustering robust standard error proposed by Chiang, Hansen, and Sasaki(2023).  The command {bf:xtregtfb} is an adaption of {bf:xtregtwo}.


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
{cmd:xtregtfb} The reported s.e. is bias-corrected CHS (BCCHS) by default. Other types of s.e. can be specified by se() where the options include chs, bchs, dka, bchs, bchsfb, dkafb. The bandwidth is chosen using Eq(6.2) and Eq(6.4) of Andrews(1991) assuming AR(1) process for {V_at} with parameters (rho_a,1-rho_a^2); the bandwidth can also be specified by lag() where any integer between 1 and the time sample size is allowed (by definition, 1 corresponds to CGM s.e.). If se() is specified with bchsfb or dkafb, the simulated fixed-b critical values (for the t-statistics under the null) will be reported at the 97.5% quantile by default. Other quantiles can be specified by level() where the default is 0.05. The critical values are regressor-specific and can be specified by whichvar() where the option is the sequence order among indepvar; the default is 1. cvsim() specifies the number of replications for simulating the fixed-b critical values; the default is 2000. bm() specifies the number of increment for approximating the standard Wiener process; the default is 1000. Other options include noconstant and fe.




{title:Reference}

{p 4 8}K. Chen and T.J. Vogelsang (2023) Fixed-b Asymptotics for Panel Models with Two-Way Clustering. Working Paper.
{p_end}

