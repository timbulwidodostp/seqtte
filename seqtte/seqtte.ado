*! version 0.11.0  03jul2026  Tom Palmer
program seqtte, eclass
    version 16

    syntax varlist(min=1 max=1 numeric) [if] [in], ///
        id(varname numeric) ///
        time(varname numeric) ///
        treatment(varname numeric) ///
        [covariates(varlist fv) ///
         ESTIMator(string) ///
         wdenominator(varlist fv) ///
         wnumerator(varlist fv) ///
         TRUNCation(real 25) ///
         SELECTIONrandom ///
         SELECTIONsample(real 0.5) ///
         SEEd(integer -1) ///
         BOOTstrap(integer 0) ///
         PLOT ///
         SURVivalmax(integer -1) ///
         expandonly]

    local outcome `varlist'

    // Default and validate estimator
    if "`estimator'" == "" local estimator "itt"
    local estimator = lower("`estimator'")
    if !inlist("`estimator'", "itt", "pp") {
        di as err "estimator() must be itt or pp"
        exit 198
    }
    // weighted_pp: censoring + IPCW weights; unweighted pp: censoring only
    local weighted_pp = ("`estimator'" == "pp" & "`wdenominator'" != "")

    // Validate selection_random option
    if "`selectionrandom'" != "" {
        if `selectionsample' <= 0 | `selectionsample' > 1 {
            di as err "selectionsample() must be in (0, 1]"
            exit 198
        }
    }

    // Validate bootstrap option
    if `bootstrap' < 0 {
        di as err "bootstrap() must be a non-negative integer"
        exit 198
    }
    local do_bs = (`bootstrap' > 0)

    // Validate survivalmax option (plot follow-up cap; -1 = adaptive default)
    if `survivalmax' != -1 & `survivalmax' < 1 {
        di as err "survivalmax() must be a positive follow-up time"
        exit 198
    }

    // expandonly returns the expanded dataset without fitting any model, so it
    // is incompatible with the analysis-dependent options.
    if "`expandonly'" != "" {
        if `do_bs' {
            di as err "expandonly cannot be combined with bootstrap()"
            exit 198
        }
        if "`plot'" != "" {
            di as err "expandonly cannot be combined with plot"
            exit 198
        }
    }

    // Capture variable count before tempvars are created
    local k_orig = c(k)
    local seed_set 0

    preserve

    marksample touse
    // markout needs base variable names, not factor-variable notation
    local cov_base `covariates'
    if "`cov_base'" != "" {
        fvrevar `cov_base', list
        local cov_base `r(varlist)'
    }
    markout `touse' `id' `time' `treatment' `cov_base'
    if "`estimator'" == "pp" {
        local wden_base `wdenominator'
        local wnum_base `wnumerator'
        if "`wden_base'" != "" {
            fvrevar `wden_base', list
            local wden_base `r(varlist)'
        }
        if "`wnum_base'" != "" {
            fvrevar `wnum_base', list
            local wnum_base `r(varlist)'
        }
        markout `touse' `wden_base' `wnum_base'
    }
    qui keep if `touse'

    sort `id' `time'

    qui count
    local n_orig = r(N)
    qui levelsof `id'
    local n_indiv = r(r)

    di as txt _n "Original dataset: " %12.0fc `n_orig' " observations, " `k_orig' " variables"

    // Event time per individual
    tempvar evttime evttime_max
    qui gen double `evttime' = .
    qui replace `evttime' = `time' if `outcome' == 1
    qui by `id': egen double `evttime_max' = max(`evttime')

    // Eligibility: not treated in any earlier period
    tempvar A_lag eligible
    qui by `id': gen byte `A_lag' = `treatment'[_n-1]
    qui gen byte `eligible' = 1
    qui replace `eligible' = 0 if `A_lag' == 1
    qui by `id': replace `eligible' = 0 ///
        if `eligible'[_n-1] == 0 & `id' == `id'[_n-1]

    qui count if `eligible' == 1
    local n_elig = r(N)
    di as txt "Eligible observations: " %12.0fc `n_elig' " (" `n_indiv' " individuals)"

    // Trial = calendar time at trial entry
    tempvar trial
    qui gen long `trial' = `time'

    // ----- PP: censoring snapshots and weight models (before expand) -----
    if "`estimator'" == "pp" {

        qui sum `time'
        local t_min = r(min)
        local t_max = r(max)
        local loop_start = `t_min' + 1

        // Per-calendar-period censoring snapshots (all PP variants)
        forvalues i = `loop_start'/`t_max' {
            tempvar _c`i'
            qui gen byte `_c`i'' = 0
            qui replace `_c`i'' = 1 ///
                if `trial' == `i' & `treatment' != `treatment'[_n-1] ///
                & `id' == `id'[_n-1]
        }
        forvalues i = `loop_start'/`t_max' {
            tempvar _ca`i'
            qui by `id': egen byte `_ca`i'' = max(`_c`i'')
        }

        // Weight models and weight snapshots (weighted PP only)
        if `weighted_pp' {
            di as txt _n "Fitting weight models..."

            // Denominator: P(A | A_lag, time + time² + time³, wdenominator)
            // Guard each fit: if treatment does not vary within the A_lag stratum
            // (e.g. monotonic, non-reversible treatment, where everyone with
            // A_lag == 1 also has treatment == 1), logistic exits with r(2000)
            // "outcome does not vary". In that case the stratum's weight
            // contribution is 1, so set the predicted probability to the constant
            // value of treatment instead of fitting a model.
            tempvar p_d0 p_d1

            qui sum `treatment' if `A_lag' == 0, meanonly
            if r(N) == 0 | r(min) == r(max) {
                qui gen double `p_d0' = cond(r(N) == 0, 1, r(mean))
            }
            else {
                qui logistic `treatment' ///
                    c.`time' c.`time'#c.`time' c.`time'#c.`time'#c.`time' ///
                    `wdenominator' if `A_lag' == 0
                qui predict double `p_d0'
            }

            qui sum `treatment' if `A_lag' == 1, meanonly
            if r(N) == 0 | r(min) == r(max) {
                qui gen double `p_d1' = cond(r(N) == 0, 1, r(mean))
            }
            else {
                qui logistic `treatment' ///
                    c.`time' c.`time'#c.`time' c.`time'#c.`time'#c.`time' ///
                    `wdenominator' if `A_lag' == 1
                qui predict double `p_d1'
            }

            // Numerator: P(A | A_lag, time + time² + time³, wnumerator)
            local stabilized = ("`wnumerator'" != "")
            if `stabilized' {
                tempvar p_n0 p_n1

                qui sum `treatment' if `A_lag' == 0, meanonly
                if r(N) == 0 | r(min) == r(max) {
                    qui gen double `p_n0' = cond(r(N) == 0, 1, r(mean))
                }
                else {
                    qui logistic `treatment' ///
                        c.`time' c.`time'#c.`time' c.`time'#c.`time'#c.`time' ///
                        `wnumerator' if `A_lag' == 0
                    qui predict double `p_n0'
                }

                qui sum `treatment' if `A_lag' == 1, meanonly
                if r(N) == 0 | r(min) == r(max) {
                    qui gen double `p_n1' = cond(r(N) == 0, 1, r(mean))
                }
                else {
                    qui logistic `treatment' ///
                        c.`time' c.`time'#c.`time' c.`time'#c.`time'#c.`time' ///
                        `wnumerator' if `A_lag' == 1
                    qui predict double `p_n1'
                }
            }

            // IPW weights (reset to 1 at each person's first record)
            tempvar ipw
            if `stabilized' {
                qui gen double `ipw' = (1 - `p_n0') / (1 - `p_d0') if `treatment' == 0
                qui replace `ipw' = `p_n1' / `p_d1'                if `treatment' == 1
            }
            else {
                qui gen double `ipw' = 1 / (1 - `p_d0') if `treatment' == 0
                qui replace `ipw' = 1 / `p_d1'          if `treatment' == 1
            }
            qui replace `ipw' = 1 if `id' != `id'[_n-1]

            // Per-calendar-period weight snapshots
            forvalues i = `loop_start'/`t_max' {
                tempvar _w`i'
                qui gen double `_w`i'' = 1
                qui replace `_w`i'' = `ipw' if `trial' == `i'
            }
            forvalues i = `loop_start'/`t_max' {
                tempvar _wa`i'
                qui by `id': egen double `_wa`i'' = max(`_w`i'')
            }

            di as txt "Weight models fitted"
        }
    }

    // ----- Expand each row to cover remaining follow-up -----
    di as txt _n "Expanding data..."

    tempvar max_t n_expand
    qui by `id': gen long `max_t' = `time'[_N]
    qui gen long `n_expand' = `max_t' - `time' + 1
    qui expand `n_expand'

    sort `id' `trial'

    tempvar time_in_trial
    qui by `id' `trial': gen long `time_in_trial' = _n

    qui drop if `eligible' == 0

    qui count
    local n_exp = r(N)
    di as txt "Expanded dataset: " %12.0fc `n_exp' " observations (person-trial periods)"

    // ----- PP: fill censoring and cumulative weights into expanded data -----
    if "`estimator'" == "pp" {
        tempvar censored
        qui gen byte `censored' = 0

        forvalues i = `loop_start'/`t_max' {
            local cond "`i' == `time_in_trial' + `trial' - 1 & `time_in_trial' > 1"
            qui replace `censored' = `_ca`i'' if `cond'
        }

        // Propagate censoring forward within each (id, trial)
        sort `id' `trial' `time_in_trial'
        qui replace `censored' = 1 ///
            if `censored'[_n-1] == 1 ///
            & `id'    == `id'[_n-1] ///
            & `trial' == `trial'[_n-1]

        if `weighted_pp' {
            // Cumulative product of weights within (id, trial)
            tempvar wt wt_cum
            qui gen double `wt' = 1
            forvalues i = `loop_start'/`t_max' {
                local cond "`i' == `time_in_trial' + `trial' - 1 & `time_in_trial' > 1"
                qui replace `wt' = `_wa`i'' if `cond'
            }
            qui gen double `wt_cum' = `wt'
            qui by `id' `trial': ///
                replace `wt_cum' = `wt_cum' * `wt_cum'[_n-1] if _n > 1

            // Truncate extreme weights
            qui replace `wt_cum' = `truncation' ///
                if `wt_cum' > `truncation' & !missing(`wt_cum')
        }

        qui count if `censored' == 0
        local n_pp = r(N)
        di as txt "Post-censoring dataset: " %12.0fc `n_pp' " observations (non-censored person-trial periods)"
    }

    // ----- Random selection of control-arm (id, trial) pairs -----
    if "`selectionrandom'" != "" {

        if `seed' != -1 & !`seed_set' {
            set seed `seed'
            local seed_set 1
        }

        // One random draw per (id, trial) pair; treatment is constant within pairs
        tempvar rnd_pair
        qui bysort `id' `trial' (`time_in_trial'): ///
            gen double `rnd_pair' = runiform() if _n == 1
        qui bysort `id' `trial' (`time_in_trial'): ///
            replace `rnd_pair' = `rnd_pair'[1]

        // Keep all treated-arm pairs; Bernoulli-sample control-arm pairs
        qui keep if `treatment' == 1 | (`treatment' == 0 & `rnd_pair' <= `selectionsample')

        qui count
        local n_sel = r(N)
        di as txt "After random selection (" %5.3f `selectionsample' ///
            " of control-arm pairs): " %12.0fc `n_sel' " observations"
    }

    // ----- Outcome and follow-up time -----
    tempvar event fu_time
    qui gen byte `event' = 0
    qui replace `event' = 1 ///
        if `evttime_max' == `time_in_trial' + `trial' - 1 ///
        & !missing(`evttime_max')
    qui gen long `fu_time' = `time_in_trial' - 1

    // ----- Readable names for the regression table -----
    // The data are preserved/restored, so these renames are discarded after the
    // fit. If a name is already in use in the caller's data we keep the
    // temporary name to avoid a collision.
    capture confirm new variable event
    if !_rc {
        rename `event' event
        local event event
    }
    capture confirm new variable followup
    if !_rc {
        rename `fu_time' followup
        local fu_time followup
    }
    capture confirm new variable trial
    if !_rc {
        rename `trial' trial
        local trial trial
    }

    // ----- expandonly: return the expanded dataset, skip the analysis -----
    if "`expandonly'" != "" {
        // period = calendar time within each trial (trial entry + follow-up)
        capture confirm new variable period
        if !_rc qui gen long period = `trial' + `fu_time'

        // Keep the per-period censoring flag (and weights) for per-protocol
        if "`estimator'" == "pp" {
            capture confirm new variable censored
            if !_rc {
                rename `censored' censored
                local censored censored
            }
            if `weighted_pp' {
                capture confirm new variable weight
                if !_rc rename `wt_cum' weight
            }
        }

        di as txt _n "expandonly: returning expanded dataset (" ///
            %12.0fc `n_exp' " observations); analysis steps skipped"

        // Post the expansion summary to e() (no model is fitted)
        ereturn clear
        ereturn scalar N_exp    = `n_exp'
        ereturn scalar N_indiv  = `n_indiv'
        ereturn scalar N_orig   = `n_orig'
        if "`selectionrandom'" != "" ereturn scalar N_sel = `n_sel'
        ereturn local  estimator "`estimator'"
        ereturn local  cmd       "seqtte"

        // restore, not cancels the automatic restore scheduled by preserve,
        // leaving the expanded sequential-trial dataset in memory. Temporary
        // variables are dropped automatically when the program ends.
        restore, not
        exit
    }

    // ----- Follow-up counts per treatment arm -----
    tempvar _fu_in0 _fu_in1 _fu_id1st
    if "`estimator'" == "pp" {
        qui bysort `id': egen byte `_fu_in0' = ///
            max(cond(`censored' == 0 & `treatment' == 0, 1, 0))
        qui bysort `id': egen byte `_fu_in1' = ///
            max(cond(`censored' == 0 & `treatment' == 1, 1, 0))
        qui count if `censored' == 0 & `treatment' == 0
        local fu_nonuniq_0 = r(N)
        qui count if `censored' == 0 & `treatment' == 1
        local fu_nonuniq_1 = r(N)
    }
    else {
        qui bysort `id': egen byte `_fu_in0' = max(`treatment' == 0)
        qui bysort `id': egen byte `_fu_in1' = max(`treatment' == 1)
        qui count if `treatment' == 0
        local fu_nonuniq_0 = r(N)
        qui count if `treatment' == 1
        local fu_nonuniq_1 = r(N)
    }
    qui bysort `id': gen byte `_fu_id1st' = (_n == 1)
    qui count if `_fu_id1st' & `_fu_in0'
    local fu_uniq_0 = r(N)
    qui count if `_fu_id1st' & `_fu_in1'
    local fu_uniq_1 = r(N)

    di as txt _n "Follow-up by treatment arm (arm 0 / arm 1):"
    di as txt "  Intervals (nonunique):  " ///
        %12.0fc `fu_nonuniq_0' " / " %12.0fc `fu_nonuniq_1'
    di as txt "  Individuals (unique):   " ///
        %12.0fc `fu_uniq_0' " / " %12.0fc `fu_uniq_1'

    local bs_ok = 0

    // ----- Cumulative-incidence follow-up cap (plot only) -----
    // Cap the plotted follow-up so the g-computation curves are not driven by
    // model extrapolation into the sparse tail. Default: the largest follow-up
    // where at least 10% of baseline trials still contribute observed
    // (uncensored) data (minimum 5 trials). Override with survivalmax().
    // Computed here, before the bootstrap, so the bands and the point curves
    // share one grid.
    if "`plot'" != "" {
        if `survivalmax' != -1 {
            local _max_fu = `survivalmax'
        }
        else {
            tempvar _cnt_ft
            if "`estimator'" == "pp" {
                qui bysort `fu_time': egen long `_cnt_ft' = total(`censored' == 0)
            }
            else {
                qui bysort `fu_time': gen long `_cnt_ft' = _N
            }
            qui sum `_cnt_ft' if `fu_time' == 0, meanonly
            local _cap_thresh = max(5, floor(r(mean) * 0.10))
            qui sum `fu_time' if `_cnt_ft' >= `_cap_thresh', meanonly
            if r(N) > 0 local _max_fu = r(max)
            else {
                qui sum `fu_time', meanonly
                local _max_fu = r(max)
            }
            qui drop `_cnt_ft'
        }
        di as txt _n "Cumulative-incidence follow-up capped at " `_max_fu'
    }

    // ----- Bootstrap -----
    if `do_bs' {

        if `seed' != -1 & !`seed_set' {
            set seed `seed'
            local seed_set 1
        }

        // Resample as many clusters as remain in the analysis sample.
        // selectionrandom (and expansion) can drop whole ids, so the original
        // individual count may exceed the clusters present here; bsample errors
        // (r(498)) if asked for more clusters than exist.
        tempvar bs_grp
        qui egen long `bs_grp' = group(`id')
        qui summarize `bs_grp', meanonly
        local n_clust = r(max)
        drop `bs_grp'

        // Save fully processed dataset (weights, outcome, renaming applied)
        tempfile bsdata
        qui save `bsdata'

        // Allocate bootstrap CIF matrices over the pre-computed follow-up grid
        // (the same cap the point curves use), so each replicate projects its
        // resampled baselines over the same 0..`_max_fu' grid.
        if "`plot'" != "" {
            tempfile _bs_grid
            local _n_t_bs = `_max_fu' + 1
            tempname bs_cif0 bs_cif1
            matrix `bs_cif0' = J(`bootstrap', `_n_t_bs', .)
            matrix `bs_cif1' = J(`bootstrap', `_n_t_bs', .)
        }

        di as txt _n "Running " `bootstrap' " bootstrap replicates..."

        tempname bs_b
        matrix `bs_b' = J(`bootstrap', 1, .)
        local bs_ok = 0
        // newid declared once; bsample creates it fresh each iteration after reload
        tempvar bs_newid

        forvalues b = 1/`bootstrap' {
            // reload from tempfile instead of nested preserve (r(621))
            qui use `bsdata', clear
            cap {
                bsample `n_clust', cluster(`id') idcluster(`bs_newid')

                if "`estimator'" == "itt" {
                    qui logistic `event' `treatment' ///
                        c.`fu_time'##c.`fu_time' ///
                        c.`trial'##c.`trial' ///
                        `covariates', cluster(`bs_newid')
                }
                else if `weighted_pp' {
                    qui logistic `event' `treatment' ///
                        c.`fu_time'##c.`fu_time' ///
                        c.`trial'##c.`trial' ///
                        `covariates' ///
                        if `censored' == 0 [pweight = `wt_cum'], cluster(`bs_newid')
                }
                else {
                    qui logistic `event' `treatment' ///
                        c.`fu_time'##c.`fu_time' ///
                        c.`trial'##c.`trial' ///
                        `covariates' ///
                        if `censored' == 0, cluster(`bs_newid')
                }

                matrix `bs_b'[`b', 1] = _b[`treatment']
                local bs_ok = `bs_ok' + 1
                if "`plot'" != "" {
                    // Refit with the treatment x follow-up interaction for the CIF
                    // bands; the bootstrap HR (bs_b) above uses the PH model.
                    if "`estimator'" == "itt" {
                        qui logistic `event' `treatment' ///
                            c.`fu_time'##c.`fu_time' ///
                            c.`trial'##c.`trial' ///
                            `covariates' c.`fu_time'#c.`treatment', cluster(`bs_newid')
                    }
                    else if `weighted_pp' {
                        qui logistic `event' `treatment' ///
                            c.`fu_time'##c.`fu_time' ///
                            c.`trial'##c.`trial' ///
                            `covariates' c.`fu_time'#c.`treatment' ///
                            if `censored' == 0 [pweight = `wt_cum'], cluster(`bs_newid')
                    }
                    else {
                        qui logistic `event' `treatment' ///
                            c.`fu_time'##c.`fu_time' ///
                            c.`trial'##c.`trial' ///
                            `covariates' c.`fu_time'#c.`treatment' ///
                            if `censored' == 0, cluster(`bs_newid')
                    }
                    // Full-grid counterfactual CIF (mirrors the main g-computation):
                    // one baseline row per (bs_newid, trial), projected over the
                    // follow-up grid, with treatment set to each arm level.
                    // bs_newid is unique per resampled cluster so duplicated
                    // clusters' trials are not merged.
                    if "`estimator'" == "pp" qui keep if `censored' == 0
                    qui keep if `fu_time' == 0
                    qui expand `_n_t_bs'
                    qui bysort `bs_newid' `trial': replace `fu_time' = _n - 1
                    qui save `_bs_grid', replace

                    // Arm 1: everyone treated
                    qui replace `treatment' = 1
                    qui predict double _bsp, pr
                    qui bysort `bs_newid' `trial' (`fu_time'): ///
                        gen double _bsls = sum(ln(1 - _bsp))
                    qui gen double _bssv = exp(_bsls)
                    collapse (mean) _sv1=_bssv, by(`fu_time')
                    sort `fu_time'
                    qui count
                    local _nb = r(N)
                    forvalues _j = 1/`_nb' {
                        local _ft = `fu_time'[`_j']
                        local _col = `_ft' + 1
                        matrix `bs_cif1'[`b', `_col'] = 1 - _sv1[`_j']
                    }

                    // Arm 0: everyone untreated
                    qui use `_bs_grid', clear
                    qui replace `treatment' = 0
                    qui predict double _bsp, pr
                    qui bysort `bs_newid' `trial' (`fu_time'): ///
                        gen double _bsls = sum(ln(1 - _bsp))
                    qui gen double _bssv = exp(_bsls)
                    collapse (mean) _sv0=_bssv, by(`fu_time')
                    sort `fu_time'
                    qui count
                    local _nb = r(N)
                    forvalues _j = 1/`_nb' {
                        local _ft = `fu_time'[`_j']
                        local _col = `_ft' + 1
                        matrix `bs_cif0'[`b', `_col'] = 1 - _sv0[`_j']
                    }
                }
            }
        }

        // Reload processed data so the main regression has the correct dataset
        qui use `bsdata', clear

        // Bootstrap SE and percentile CI (95%) from log-HR distribution
        // svmat places values in obs 1..B; remaining obs are missing — both
        // sum and _pctile skip missing values automatically.
        qui svmat double `bs_b', names(_seqtte_bs_)
        qui sum _seqtte_bs_1
        local bs_se = r(sd)
        qui _pctile _seqtte_bs_1, p(2.5 97.5)
        local bs_ll = r(r1)
        local bs_ul = r(r2)
        drop _seqtte_bs_1
    }

    // ----- Pooled logistic regression -----
    di as txt _n "Fitting " upper("`estimator'") " model..."

    // Fit with logit (quietly) and display the exponentiated estimates labelled
    // as hazard ratios. The pooled logistic regression approximates a discrete-
    // time hazard model, so exp(coef) is reported as a hazard ratio to match the
    // R and Python implementations of this method.
    if "`estimator'" == "itt" {
        qui logit `event' `treatment' ///
            c.`fu_time'##c.`fu_time' ///
            c.`trial'##c.`trial' ///
            `covariates', cluster(`id')
    }
    else if `weighted_pp' {
        qui logit `event' `treatment' ///
            c.`fu_time'##c.`fu_time' ///
            c.`trial'##c.`trial' ///
            `covariates' ///
            if `censored' == 0 [pweight = `wt_cum'], cluster(`id')
    }
    else {
        qui logit `event' `treatment' ///
            c.`fu_time'##c.`fu_time' ///
            c.`trial'##c.`trial' ///
            `covariates' ///
            if `censored' == 0, cluster(`id')
    }

    _coef_table_header
    ereturn display, eform("Haz. ratio")

    // Bootstrap summary (after the fit so the point HR is available)
    if `do_bs' {
        local bs_hr    = exp(_b[`treatment'])
        local bs_hr_ll = exp(`bs_ll')
        local bs_hr_ul = exp(`bs_ul')
        di as txt _n "Bootstrap complete: " `bs_ok' "/" `bootstrap' " replicates succeeded"
        di as txt "Bootstrap SE (log-HR):      " %7.4f `bs_se'
        di as txt "Bootstrap 95% CI (log-HR): [" %7.4f `bs_ll' ", " %7.4f `bs_ul' "]"
        local bs_ci_ll = ltrim(string(`bs_hr_ll', "%7.4f"))
        di as txt "Estimated HR: " %7.4f `bs_hr' ///
            " with bootstrap 95% CI [" "`bs_ci_ll'" ", " %7.4f `bs_hr_ul' "]"
    }

    // ----- Cumulative incidence by g-computation -----
    // G-formula: for each person-trial compute the individual-level cumulative
    // survival S_i(t) = prod_{s<=t}(1 - h_i(s)), then average S_i(t) across
    // all person-trials at each follow-up time. CIF = 1 - mean(S_i(t)).
    // This matches the R SEQTaRget implementation and avoids the Jensen's
    // inequality distortion of the "average hazard then product-limit" approach.
    // For weighted PP, IPCW weights are applied when averaging arm-0 survival
    // to correct for informative censoring.
    // Follow-up is truncated where arm-1 has fewer than 5 person-trials to
    // avoid unstable estimates in the sparse tail.
    if "`plot'" != "" {
        tempvar _pred _logsurv _surv
        // Refit the outcome model with a treatment x follow-up interaction so the
        // cumulative-incidence curves let the treatment effect vary over follow-up
        // (matching R SEQTaRget / Python pySEQTarget, which add haart_bas*followup
        // for the survival curves). The reported hazard ratio keeps the
        // proportional-hazards (no-interaction) model fitted above, so we store it
        // and restore it after the g-computation.
        tempname _est_ph
        qui estimates store `_est_ph'
        if "`estimator'" == "itt" {
            qui logit `event' `treatment' ///
                c.`fu_time'##c.`fu_time' ///
                c.`trial'##c.`trial' ///
                `covariates' c.`fu_time'#c.`treatment', cluster(`id')
        }
        else if `weighted_pp' {
            qui logit `event' `treatment' ///
                c.`fu_time'##c.`fu_time' ///
                c.`trial'##c.`trial' ///
                `covariates' c.`fu_time'#c.`treatment' ///
                if `censored' == 0 [pweight = `wt_cum'], cluster(`id')
        }
        else {
            qui logit `event' `treatment' ///
                c.`fu_time'##c.`fu_time' ///
                c.`trial'##c.`trial' ///
                `covariates' c.`fu_time'#c.`treatment' ///
                if `censored' == 0, cluster(`id')
        }
        // Counterfactual g-computation, standardised to the baseline covariate
        // distribution of ALL person-trials. Keeping the observed treated /
        // control subsets ("keep if treatment==1/0") standardises each arm to a
        // different, confounded covariate distribution; instead we set treatment
        // to each arm level for every person-trial and predict. Each trial's
        // baseline (fu_time==0) covariates are held constant (the expand step
        // duplicated the trial-entry row), so we project each baseline over the
        // full follow-up grid, predict the hazard from the interacted (IPCW-
        // weighted) model, form S_i(t) = prod_{s<=t}(1 - h_i(s)), then average
        // S_i(t) across all trials at each t with a simple mean (the IPCW is
        // already in the model fit). CIF = 1 - mean(S_i(t)). Matches R/Python.
        tempfile _cif_base _grid _arm1_data
        qui save `_cif_base'

        // Follow-up grid length from the pre-computed cap (survivalmax / adaptive).
        local _n_t = `_max_fu' + 1
        di as txt _n "Cumulative incidence by g-computation over follow-up 0-" `_max_fu'

        // One baseline (fu_time==0, never censored) row per (id, trial),
        // replicated over the full follow-up grid.
        qui keep if `fu_time' == 0
        qui expand `_n_t'
        qui bysort `id' `trial': replace `fu_time' = _n - 1
        qui save `_grid'

        // --- Arm 1: everyone treated ---
        qui replace `treatment' = 1
        qui predict double `_pred', pr
        qui bysort `id' `trial' (`fu_time'): ///
            gen double `_logsurv' = sum(ln(1 - `_pred'))
        qui gen double `_surv' = exp(`_logsurv')
        collapse (mean) _msurv1=`_surv', by(`fu_time')
        sort `fu_time'
        qui gen double _cif1 = 1 - _msurv1
        qui save `_arm1_data'

        // --- Arm 0: everyone untreated ---
        qui use `_grid', clear
        qui replace `treatment' = 0
        qui predict double `_pred', pr
        qui bysort `id' `trial' (`fu_time'): ///
            gen double `_logsurv' = sum(ln(1 - `_pred'))
        qui gen double `_surv' = exp(`_logsurv')
        collapse (mean) _msurv0=`_surv', by(`fu_time')
        sort `fu_time'
        qui gen double _cif0 = 1 - _msurv0

        // --- Combine ---
        qui merge 1:1 `fu_time' using `_arm1_data', nogen
        sort `fu_time'
        qui count
        local _n_t = r(N)
        tempname _cif_mat
        matrix `_cif_mat' = J(`_n_t', 3, .)
        matrix colnames `_cif_mat' = fu_time cif0 cif1
        forvalues _i = 1/`_n_t' {
            matrix `_cif_mat'[`_i', 1] = `fu_time'[`_i']
            matrix `_cif_mat'[`_i', 2] = _cif0[`_i']
            matrix `_cif_mat'[`_i', 3] = _cif1[`_i']
        }

        // --- Bootstrap CI bands ---
        if `do_bs' & `bs_ok' > 0 {
            local _n_t_bs = `_max_fu' + 1
            clear
            qui set obs `bootstrap'
            qui svmat double `bs_cif0', names(_c0_)
            qui svmat double `bs_cif1', names(_c1_)
            tempname _ci0_lo _ci0_hi _ci1_lo _ci1_hi
            matrix `_ci0_lo' = J(`_n_t_bs', 1, .)
            matrix `_ci0_hi' = J(`_n_t_bs', 1, .)
            matrix `_ci1_lo' = J(`_n_t_bs', 1, .)
            matrix `_ci1_hi' = J(`_n_t_bs', 1, .)
            forvalues _j = 1/`_n_t_bs' {
                capture qui _pctile _c0_`_j', p(2.5 97.5)
                if _rc == 0 {
                    matrix `_ci0_lo'[`_j', 1] = r(r1)
                    matrix `_ci0_hi'[`_j', 1] = r(r2)
                }
                capture qui _pctile _c1_`_j', p(2.5 97.5)
                if _rc == 0 {
                    matrix `_ci1_lo'[`_j', 1] = r(r1)
                    matrix `_ci1_hi'[`_j', 1] = r(r2)
                }
            }
            tempname _cif_mat_ci
            matrix `_cif_mat_ci' = J(`_n_t', 7, .)
            matrix colnames `_cif_mat_ci' = fu_time cif0 cif0_lo cif0_hi cif1 cif1_lo cif1_hi
            forvalues _i = 1/`_n_t' {
                local _ft  = round(`_cif_mat'[`_i', 1])
                local _col = `_ft' + 1
                matrix `_cif_mat_ci'[`_i', 1] = `_cif_mat'[`_i', 1]
                matrix `_cif_mat_ci'[`_i', 2] = `_cif_mat'[`_i', 2]
                matrix `_cif_mat_ci'[`_i', 3] = `_ci0_lo'[`_col', 1]
                matrix `_cif_mat_ci'[`_i', 4] = `_ci0_hi'[`_col', 1]
                matrix `_cif_mat_ci'[`_i', 5] = `_cif_mat'[`_i', 3]
                matrix `_cif_mat_ci'[`_i', 6] = `_ci1_lo'[`_col', 1]
                matrix `_cif_mat_ci'[`_i', 7] = `_ci1_hi'[`_col', 1]
            }
            matrix `_cif_mat' = `_cif_mat_ci'
        }

        qui use `_cif_base', clear
        // Restore the proportional-hazards model so the posted e(b)/e(V) and the
        // reported hazard ratio are unaffected by the interacted plot model.
        qui estimates restore `_est_ph'
        qui estimates drop `_est_ph'
    }

    restore

    ereturn scalar N_indiv   = `n_indiv'
    ereturn scalar N_orig    = `n_orig'
    ereturn scalar N_exp     = `n_exp'
    if "`selectionrandom'" != "" {
        ereturn scalar N_sel            = `n_sel'
        ereturn scalar selection_sample = `selectionsample'
    }
    if `do_bs' {
        ereturn scalar N_boot = `bs_ok'
        ereturn scalar bs_se  = `bs_se'
        ereturn scalar bs_ll  = `bs_ll'
        ereturn scalar bs_ul  = `bs_ul'
        ereturn matrix bs_b   = `bs_b'
    }
    ereturn scalar N_nonuniq_arm0 = `fu_nonuniq_0'
    ereturn scalar N_nonuniq_arm1 = `fu_nonuniq_1'
    ereturn scalar N_uniq_arm0    = `fu_uniq_0'
    ereturn scalar N_uniq_arm1    = `fu_uniq_1'
    if "`plot'" != "" {
        ereturn matrix cif = `_cif_mat'
    }
    ereturn local  estimator "`estimator'"
    ereturn local  cmd       "seqtte"

    // ----- Cumulative incidence plot -----
    if "`plot'" != "" {
        preserve
        clear
        tempname _tmp
        matrix `_tmp' = e(cif)
        local _nr = rowsof(`_tmp')
        local _nc = colsof(`_tmp')
        qui set obs `_nr'
        qui svmat double `_tmp', names(col)
        // Arm colours are taken from the current graph scheme rather than
        // hardcoded, so that monochrome schemes such as sj (see help scheme)
        // give a figure suitable for black-and-white print.  pstyle() ties
        // each confidence band to the line of the same arm, which would
        // otherwise pick up the third and fourth plot styles.  The two line
        // patterns are set explicitly so that the arms remain distinguishable
        // under a monochrome scheme, where both colours are black.
        if `_nc' == 7 {
            twoway ///
                (rarea cif0_lo cif0_hi fu_time, pstyle(p1) fintensity(20) lwidth(none)) ///
                (rarea cif1_lo cif1_hi fu_time, pstyle(p2) fintensity(20) lwidth(none)) ///
                (line cif0 fu_time, pstyle(p1) lpattern(solid) lwidth(medthick)) ///
                (line cif1 fu_time, pstyle(p2) lpattern(dash)  lwidth(medthick)), ///
                ytitle("Cumulative incidence") ///
                xtitle("Follow-up time") ///
                title("Cumulative incidence by treatment arm") ///
                legend(order(3 "Arm 0 (control)" 4 "Arm 1 (treated)")) ///
                name(seqtte_cif, replace)
        }
        else {
            twoway ///
                (line cif0 fu_time, pstyle(p1) lpattern(solid) lwidth(medthick)) ///
                (line cif1 fu_time, pstyle(p2) lpattern(dash)  lwidth(medthick)), ///
                ytitle("Cumulative incidence") ///
                xtitle("Follow-up time") ///
                title("Cumulative incidence by treatment arm") ///
                legend(label(1 "Arm 0 (control)") label(2 "Arm 1 (treated)")) ///
                name(seqtte_cif, replace)
        }
        restore
    }
end
