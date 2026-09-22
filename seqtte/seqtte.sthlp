{smcl}
{* *! version 0.11.0  03jul2026  Tom Palmer}{...}
{vieweralsosee "seqtte" "help seqtte"}{...}
{viewerjumpto "Syntax" "seqtte##syntax"}{...}
{viewerjumpto "Description" "seqtte##description"}{...}
{viewerjumpto "Options" "seqtte##options"}{...}
{viewerjumpto "Examples" "seqtte##examples"}{...}
{viewerjumpto "Stored results" "seqtte##results"}{...}
{viewerjumpto "References" "seqtte##references"}{...}
{viewerjumpto "Author" "seqtte##author"}{...}
{title:Title}

{phang}
{bf:seqtte} {hline 2} Sequential target trial emulation

{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:seqtte} {it:outcomevar} {ifin}{cmd:,}
{cmdab:id(}{it:varname}{cmd:)}
{cmdab:time(}{it:varname}{cmd:)}
{cmdab:treatment(}{it:varname}{cmd:)}
[{it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt id(varname)}}individual identifier variable{p_end}
{synopt:{opt time(varname)}}integer calendar time variable{p_end}
{synopt:{opt treatment(varname)}}binary treatment indicator (0/1){p_end}
{synoptline}
{syntab:Optional}
{synopt:{opt covariates(varlist)}}adjustment covariates for the outcome model{p_end}
{synopt:{opt estimator(string)}}{cmd:itt} (default) or {cmd:pp}{p_end}
{synopt:{opt wdenominator(varlist)}}denominator weight model covariates; if omitted with {cmd:pp}, an unweighted per-protocol analysis is performed{p_end}
{synopt:{opt wnumerator(varlist)}}numerator weight model covariates; if supplied, stabilized weights are used{p_end}
{synopt:{opt truncation(#)}}truncation threshold for cumulative weights; default 25{p_end}
{synopt:{opt selectionrandom}}randomly subsample the control-arm (id, trial) pairs{p_end}
{synopt:{opt selectionsample(#)}}proportion of control-arm pairs to retain under {cmd:selectionrandom}; default 0.5{p_end}
{synopt:{opt bootstrap(#)}}number of bootstrap replicates for the standard error and 95% percentile CI; default 0{p_end}
{synopt:{opt seed(#)}}random-number seed; default {cmd:-1} (seed not set){p_end}
{synopt:{opt plot}}plot cumulative incidence curves for each treatment arm{p_end}
{synopt:{opt survivalmax(#)}}cap the follow-up shown in the cumulative incidence curves; default is adaptive{p_end}
{synopt:{opt expandonly}}return the expanded sequential-trial dataset and skip the analysis{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:seqtte} estimates the causal effect of a sustained treatment strategy
using sequential target trial emulation
({help seqtte##HR2016:Hernán and Robins, 2016}).
Two estimators are available.

{pstd}
{ul:Intent-to-treat (ITT)} ({cmd:estimator(itt)}, the default).
Estimates the effect of being assigned to treatment at trial entry,
regardless of subsequent treatment changes.
No weighting is used.

{pstd}
{ul:Per-protocol (PP)} ({cmd:estimator(pp)}).
Estimates the effect of sustained adherence to the assigned treatment strategy.
Individuals are censored at the period in which they deviate from their
assigned treatment.
If {cmd:wdenominator()} is supplied, inverse probability of censoring weights
(IPCW) are applied to adjust for informative censoring (weighted PP).
If {cmd:wdenominator()} is omitted, censoring is applied but no weights are
used (unweighted PP).

{pstd}
{ul:Input data.}
The data should be in long (person-period) format with one row per individual
per time period.
{it:outcomevar} should equal 1 only in the period the event first occurs and
0 otherwise.
The {it:time} variable should take consecutive integer values.

{pstd}
{ul:Algorithm.}
For each eligible person-period (i.e. periods in which the individual has
not yet received treatment), a trial is initiated.
The individual is then followed from that trial entry time to the end of
their observed follow-up.
Data are expanded so that each person can contribute to multiple trials.
A pooled logistic regression model is then fitted with quadratic polynomial
terms for follow-up time within trial and trial number, with standard errors
clustered by individual.
The pooled logistic regression approximates a discrete-time hazard model, so
the exponentiated treatment coefficient is reported as a hazard ratio.

{pstd}
{ul:Weight models (weighted PP only).}
When {cmd:wdenominator()} is supplied with {cmd:estimator(pp)},
four logistic regression models are fitted on the pre-expansion data,
stratified by prior treatment status ({it:A_lag} = 0 or 1):
a denominator model including {cmd:wdenominator()} covariates
and a cubic polynomial in calendar time,
and (if {cmd:wnumerator()} is supplied) a numerator model including
{cmd:wnumerator()} covariates and the same time polynomial.
Unstabilized weights are used when {cmd:wnumerator()} is omitted;
stabilized weights (numerator/denominator) are used when it is supplied.
Cumulative products of per-period weights are formed within each trial,
then truncated at {cmd:truncation()}.
When {cmd:wdenominator()} is omitted, no weight models are fitted and the
outcome model is fitted on the censored data without weighting.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt id(varname)} specifies the variable identifying individuals.

{phang}
{opt time(varname)} specifies the integer calendar time variable.
Consecutive integer values represent consecutive periods.

{phang}
{opt treatment(varname)} specifies the binary treatment indicator
(0 = untreated, 1 = treated).

{dlgtab:Optional}

{phang}
{opt covariates(varlist)} specifies adjustment covariates for the
outcome model.
In the expanded dataset each person-trial record takes covariate
values from the trial entry period, so both time-fixed and
time-varying covariates are included at their trial-entry values.
{help fvvarlist:Factor-variable notation} (e.g. {cmd:i.group}) is allowed.

{phang}
{opt estimator(string)} specifies the estimator: {cmd:itt} (default)
for the intent-to-treat effect, or {cmd:pp} for the per-protocol effect.

{phang}
{opt wdenominator(varlist)} specifies the covariates for the denominator
weight models.
These are fitted on the pre-expansion data and should include all
time-varying confounders of the treatment–outcome relationship.
When supplied with {cmd:estimator(pp)}, IPCW weights are applied (weighted PP);
when omitted, censoring is applied without weighting (unweighted PP).
{help fvvarlist:Factor-variable notation} (e.g. {cmd:i.group}) is allowed.

{phang}
{opt wnumerator(varlist)} specifies the covariates for the numerator
(stabilization) weight models.
These should be a subset of {cmd:wdenominator()}, typically restricted to
baseline (study-entry) values of covariates.
When omitted, unstabilized weights are used.
{help fvvarlist:Factor-variable notation} (e.g. {cmd:i.group}) is allowed.

{phang}
{opt truncation(#)} specifies the upper truncation threshold applied
to the cumulative IPW weights.
Default is 25.

{phang}
{opt selectionrandom} randomly subsamples the control-arm (id, trial) pairs
(Bernoulli sampling), retaining all treated-arm pairs.
This reduces the size of the expanded dataset and is useful when the expanded
data are very large.
By default no subsampling is applied (all pairs are used).

{phang}
{opt selectionsample(#)} specifies the proportion of control-arm pairs to
retain when {cmd:selectionrandom} is specified.
Must be in (0, 1]; default is 0.5.
It has no effect unless {cmd:selectionrandom} is also specified.

{phang}
{opt bootstrap(#)} requests {it:#} bootstrap replicates (resampling individuals
with replacement) to compute a bootstrap standard error and 95% percentile
confidence interval for the log-hazard ratio, reported below the regression
table.
Default is 0 (no bootstrap).

{phang}
{opt seed(#)} sets the random-number seed, for reproducibility of the bootstrap
resampling and of {cmd:selectionrandom}.
Default is {cmd:-1}, meaning the seed is not set.

{phang}
{opt plot} produces cumulative incidence (1 - survival) curves for each
treatment arm by g-computation from the fitted outcome model.
Cannot be combined with {cmd:expandonly}.

{phang}
{opt survivalmax(#)} caps the follow-up time shown in the cumulative incidence
curves (with {cmd:plot}).
The g-computation projects every trial over the follow-up grid, so at long
follow-up times, where few trials contributed observed data, the curves are
driven by model extrapolation.
By default the follow-up is capped at the largest time where at least 10% of
the baseline trials still contribute observed (uncensored) data (minimum 5
trials); specify {opt survivalmax(#)} to set the cap explicitly.

{phang}
{opt expandonly} performs the data expansion only and leaves the expanded
sequential-trial dataset in memory, skipping the weight models, outcome model,
bootstrap, and cumulative-incidence steps.
The returned data contain the original variables together with {cmd:trial}
(calendar time of trial entry), {cmd:followup} (time since trial entry),
{cmd:period} (calendar time, equal to {cmd:trial} + {cmd:followup}), and
{cmd:event} (the period-specific outcome indicator); for {cmd:estimator(pp)} a
{cmd:censored} indicator (and, for weighted PP, the cumulative {cmd:weight}) are
also included.
This option cannot be combined with {cmd:bootstrap()} or {cmd:plot}.

{marker examples}{...}
{title:Examples}

{pstd}Setup: generate a synthetic person-period dataset with a baseline
covariate ({cmd:age}) and a time-varying confounder ({cmd:bmi}){p_end}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. set obs 500}{p_end}
{phang2}{cmd:. gen id = _n}{p_end}
{phang2}{cmd:. gen age = rnormal(50, 10)}{p_end}
{phang2}{cmd:. expand 12}{p_end}
{phang2}{cmd:. bysort id (age): gen time = _n - 1}{p_end}
{phang2}{cmd:. gen bmi = .}{p_end}
{phang2}{cmd:. bysort id (time): replace bmi = rnormal(25, 4) if time == 0}{p_end}
{phang2}{cmd:. bysort id (time): replace bmi = bmi[_n-1] + rnormal(0, 1) if time > 0}{p_end}
{phang2}{cmd:. gen treatment = .}{p_end}
{phang2}{cmd:. bysort id (time): replace treatment = (runiform() < invlogit(-2 + 0.1 * (bmi - 25))) if time == 0}{p_end}
{phang2}{cmd:. bysort id (time): replace treatment = cond(treatment[_n-1] == 0, runiform() < invlogit(-2 + 0.1 * (bmi - 25)), runiform() < 0.7) if time > 0}{p_end}
{phang2}{cmd:. gen outcome = (runiform() < invlogit(-3 + 0.05 * (bmi - 25) - 0.4 * treatment))}{p_end}
{phang2}{cmd:. bysort id (time): gen cumev = sum(outcome)}{p_end}
{phang2}{cmd:. drop if cumev > 1}{p_end}
{phang2}{cmd:. drop cumev}{p_end}

{pstd}ITT estimator{p_end}

{phang2}{cmd:. seqtte outcome, id(id) time(time) treatment(treatment) covariates(age)}{p_end}

{pstd}Unweighted PP estimator (censoring applied, no IPCW weights){p_end}

{phang2}{cmd:. seqtte outcome, id(id) time(time) treatment(treatment) covariates(age) estimator(pp)}{p_end}

{pstd}PP estimator with unstabilized weights{p_end}

{phang2}{cmd:. seqtte outcome, id(id) time(time) treatment(treatment) covariates(age) estimator(pp) wdenominator(age bmi)}{p_end}

{pstd}PP estimator with stabilized weights. Note that the numerator model should
be a submodel of the denominator model; if the same covariates are given to both
the models coincide, every weight is 1, and the fit is identical to the
unweighted per-protocol analysis.{p_end}

{phang2}{cmd:. seqtte outcome, id(id) time(time) treatment(treatment) covariates(age) estimator(pp) wdenominator(age bmi) wnumerator(age)}{p_end}

{pstd}Cumulative incidence curves by treatment arm{p_end}

{phang2}{cmd:. seqtte outcome, id(id) time(time) treatment(treatment) covariates(age) plot}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:seqtte} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations in the pooled logistic regression{p_end}
{synopt:{cmd:e(N_indiv)}}number of individuals in the original data{p_end}
{synopt:{cmd:e(N_orig)}}number of observations in the original data{p_end}
{synopt:{cmd:e(N_exp)}}number of observations in the expanded dataset{p_end}
{synopt:{cmd:e(r2_p)}}pseudo R-squared{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(N_uniq_arm0)}}number of individuals contributing follow-up to arm 0{p_end}
{synopt:{cmd:e(N_uniq_arm1)}}number of individuals contributing follow-up to arm 1{p_end}
{synopt:{cmd:e(N_nonuniq_arm0)}}number of follow-up intervals in arm 0{p_end}
{synopt:{cmd:e(N_nonuniq_arm1)}}number of follow-up intervals in arm 1{p_end}
{synopt:{cmd:e(N_sel)}}observations after random selection (if {cmd:selectionrandom}){p_end}
{synopt:{cmd:e(selection_sample)}}control-arm sampling proportion (if {cmd:selectionrandom}){p_end}
{synopt:{cmd:e(N_boot)}}number of successful bootstrap replicates (if {cmd:bootstrap()}){p_end}
{synopt:{cmd:e(bs_se)}}bootstrap standard error of the log-hazard ratio (if {cmd:bootstrap()}){p_end}
{synopt:{cmd:e(bs_ll)}}lower bound of the bootstrap 95% percentile CI for the hazard ratio (if {cmd:bootstrap()}){p_end}
{synopt:{cmd:e(bs_ul)}}upper bound of the bootstrap 95% percentile CI for the hazard ratio (if {cmd:bootstrap()}){p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:seqtte}{p_end}
{synopt:{cmd:e(estimator)}}{cmd:itt} or {cmd:pp}{p_end}
{synopt:{cmd:e(depvar)}}name of the outcome variable{p_end}
{synopt:{cmd:e(clustvar)}}name of the cluster variable{p_end}
{synopt:{cmd:e(vcetype)}}{cmd:Robust}{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector (log scale; exponentiated values reported as hazard ratios){p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}
{synopt:{cmd:e(bs_b)}}bootstrap log-hazard-ratio replicates (if {cmd:bootstrap()}){p_end}
{synopt:{cmd:e(cif)}}cumulative incidence by arm and follow-up time (if {cmd:plot}){p_end}

{marker references}{...}
{title:References}

{marker HR2016}{...}
{phang}
Hernán MA, Robins JM. 2016.
Using Big Data to Emulate a Target Trial When a Randomized Trial Is Not Available.
{it:American Journal of Epidemiology} 183(8): 758–764.

{phang}
Danaei G, Rodríguez LAG, Cantero OF, Logan R, Hernán MA. 2013.
Observational data for comparative effectiveness research:
an emulation of randomised trials of statins and primary prevention of coronary heart disease.
{it:Statistical Methods in Medical Research} 22(1): 70–96.

{phang}
Maringe C, Benitez Majano S, Exarchakou A, et al. 2020.
Reflection on modern methods: trial emulation in the presence of immortal-time bias.
Assessing the benefit of major surgery for elderly lung cancer patients using
observational data.
{it:International Journal of Epidemiology} 49(5): 1719–1729.

{marker author}{...}
{title:Author}

{phang}
Tom Palmer, University of Bristol, Bristol, UK.
{browse "mailto:remlapmot@hotmail.com":remlapmot@hotmail.com}

{phang}
Michalis Katsoulis, UCL, London, UK.

{phang}
Please report any bugs or feature requests at
{browse "https://github.com/remlapmot/seqtte/issues"}.
{p_end}
