********************************************************************************
* Threshold Policy Learning via Strata Means                                   *
* Percentile path (IPW value) + stratified bootstrap bands + diagnostics       *
*                                                                              *
* Date: 22/09/2025                                                             *
* Author: Roberto Vacante                                                      *
********************************************************************************
version 16.0

* --------- Paths -----------------------------------------
* DATA  : folder that contains data
* OUTDIR: where figures and logs will be written
local DATA   "data"
local OUTDIR "output/TPL_results"

cap mkdir "output"
cap mkdir "`OUTDIR'"

* --------- Optional aesthetics ------------------------------------------------
cap set scheme s2color
graph set window fontface "Arial Narrow"
graph set print  fontface "Arial Narrow"

* --------- Variables ----------------------------------------------------------
local treatvar vti_treated
local bin_outcomes  m_has_job m_jobactive m_start_business m_has_savings
local cont_outcomes m_risk_general m_risk_job m_tot_savings m_monthly_income m_life_satis_now

* --------- Log ----------------------------------------------------------------
cap log close _all
local today "$S_DATE"
log using "`OUTDIR'/TPL_`today'", replace

* Export size for PNG
local EX_W 4000
local EX_H 4800

* ==============================================================================
* Main loop (binary + continuous outcomes)
* ==============================================================================
foreach outcomevar in `bin_outcomes' `cont_outcomes' {

    * Load panel and basic filters (adapt if needed)
    use "`DATA'/base&midline.dta", clear
    keep if ACAV_force == 0 & new_sampled == 0

    * Strata (female × location)
    egen strata = group(female location), label

    * Minimal set
    keep id `treatvar' `outcomevar' strata
    drop if missing(`treatvar', `outcomevar', strata)

    * Need variation in treatment overall
    quietly count if `treatvar'==1
    if r(N)==0 | r(N)==_N continue

    * Block propensities (design-known)
    bysort strata: egen double ps = mean(`treatvar')

    * Score: strata mean contrast tau_hat = mu1 - mu0
    tempvar y1only y0only y1bar y0bar
    gen double `y1only' = cond(`treatvar'==1, `outcomevar', .)
    gen double `y0only' = cond(`treatvar'==0, `outcomevar', .)
    bysort strata: egen double `y1bar' = mean(`y1only')
    bysort strata: egen double `y0bar' = mean(`y0only')
    gen double tau_hat = `y1bar' - `y0bar'
    drop if missing(tau_hat)

    * RCT mean
    quietly summarize `outcomevar'
    local RCT_mean = r(mean)

    * ---------------- Percentile path (value + diagnostics) -------------------
    tempfile tmpplot
    capture postclose __P
    postfile __P int p double W_hat RCT_mean Diff_TR_RCT ///
        long N11 long N10 long N01 long N00 ///
        double mismatch coverage using "`tmpplot'", replace

    tempvar Wi
    gen double `Wi' = .

    forvalues p = 5(5)95 {
        quietly _pctile tau_hat, p(`p')
        local t_thresh = r(r1)

        gen byte d_theta = (tau_hat > `t_thresh') if tau_hat < .

        replace `Wi' = .
        replace `Wi' = `outcomevar' * ( ///
            (`treatvar'==1 & d_theta==1)/ps + ///
            (`treatvar'==0 & d_theta==0)/(1-ps) ) if ///
            d_theta<. & ps>0 & ps<1

        quietly summarize `Wi'
        local W_hat = r(mean)
        local diff  = `W_hat' - `RCT_mean'

        quietly count if d_theta==1 & `treatvar'==1
        local N11 = r(N)
        quietly count if d_theta==1 & `treatvar'==0
        local N10 = r(N)
        quietly count if d_theta==0 & `treatvar'==1
        local N01 = r(N)
        quietly count if d_theta==0 & `treatvar'==0
        local N00 = r(N)

        local Ntot     = `N11' + `N10' + `N01' + `N00'
        local mismatch = cond(`Ntot'>0, (`N10' + `N01') / `Ntot', .)
        local coverage = cond(`Ntot'>0, (`N11' + `N10') / `Ntot', .)

        post __P (`p') (`W_hat') (`RCT_mean') (`diff') ///
                  (`N11') (`N10') (`N01') (`N00') ///
                  (`mismatch') (`coverage')

        drop d_theta
    }
    postclose __P

    * ---------------- Stratified bootstrap (bands for Diff_TR_RCT) -----------
    local reps = 500
    tempfile basepanel
    preserve
        keep id strata `treatvar' `outcomevar'
        drop if missing(strata, `treatvar', `outcomevar')
        save `basepanel', replace
    restore

    tempfile bootrows
    capture postclose __B
    postfile __B int p int rep double diff using "`bootrows'", replace

    set seed 13579
    forvalues b = 1/`reps' {
        use `basepanel', clear
        bsample, strata(strata)

        bysort strata: egen double ps = mean(`treatvar')

        tempvar y1only y0only y1bar y0bar WiB tau_b
        gen double `y1only' = cond(`treatvar'==1, `outcomevar', .)
        gen double `y0only' = cond(`treatvar'==0, `outcomevar', .)
        bysort strata: egen double `y1bar' = mean(`y1only')
        bysort strata: egen double `y0bar' = mean(`y0only')
        gen double `tau_b' = `y1bar' - `y0bar'
        drop if missing(`tau_b')

        quietly summarize `outcomevar'
        local RCTm = r(mean)

        forvalues p = 5(5)95 {
            quietly _pctile `tau_b', p(`p')
            local t_thresh = r(r1)
            gen byte d_theta = (`tau_b' > `t_thresh') if `tau_b' < .

            gen double `WiB' = .
            replace `WiB' = `outcomevar' * ( ///
                (`treatvar'==1 & d_theta==1)/ps + ///
                (`treatvar'==0 & d_theta==0)/(1-ps) ) if ///
                d_theta<. & ps>0 & ps<1

            quietly summarize `WiB'
            local W_hat = r(mean)
            local diff  = `W_hat' - `RCTm'

            post __B (`p') (`b') (`diff')
            drop d_theta `WiB'
        }
    }
    postclose __B

    * Collapse bootstrap to CI and merge with path
    use "`bootrows'", clear
    bysort p: egen diff_lo   = pctile(diff), p(2.5)
    bysort p: egen diff_hi   = pctile(diff), p(97.5)
    bysort p: egen diff_mean = mean(diff)
    keep p diff_mean diff_lo diff_hi
    duplicates drop
    tempfile bootsum
    save `bootsum', replace

    use "`tmpplot'", clear
    merge 1:1 p using `bootsum', nogenerate

    * ---------------- Plots: composition+mismatch (top) | TR−RCT (bottom) ----
    count if !missing(Diff_TR_RCT, p)
    if r(N) > 0 {

        * Scale for binary outcomes (percentage points)
        local scale = 1
        foreach b of local bin_outcomes {
            if ("`outcomevar'"=="`b'") local scale = 100
        }

        * Top panel: stacked cells + mismatch (y2)
        gen long tot = N11 + N10 + N01 + N00
        gen long b11 = 0
        gen long t11 = N11
        gen long b00 = t11
        gen long t00 = b00 + N00
        gen long b10 = t00
        gen long t10 = b10 + N10
        gen long b01 = t10
        gen long t01 = b01 + N01

        twoway ///
            (rbar b11 t11 p, barwidth(3) fcolor(navy%75)   lcolor(navy%90)) ///
            (rbar b00 t00 p, barwidth(3) fcolor(gs12%75)   lcolor(gs8%90)) ///
            (rbar b10 t10 p, barwidth(3) fcolor(maroon%30) lcolor(maroon%20)) ///
            (rbar b01 t01 p, barwidth(3) fcolor(orange%30) lcolor(orange%20)) ///
            (line mismatch p, yaxis(2) lpattern(solid) lcolor(black)) ///
            , ///
              ytitle("Observations", size(medium)) ///
              ytitle("Mismatch share", axis(2) size(medium)) ///
              ylabel(, angle(h)) ///
              ylabel(0(.25)1, axis(2) angle(h)) ///
              yscale(range(0 1) axis(2)) ///
              xtitle("Percentile threshold of score") ///
              legend(order(1 2 3 4) rows(2) cols(2) ///
                     label(1 "TR=1 & T=1") ///
                     label(2 "TR=0 & T=0") ///
                     label(3 "TR=1 & T=0") ///
                     label(4 "TR=0 & T=1")) ///
              note("Mismatch: share where TR and experimental T disagree.", size(small)) ///
              title("Group composition — `outcomevar'", size(medlarge)) ///
              graphregion(color(white)) ///
              name(G_top_`outcomevar', replace) nodraw

        * Bottom panel: TR − RCT with bootstrap band
        gen double diff_s    = Diff_TR_RCT * `scale'
        gen double diff_lo_s = diff_lo     * `scale'
        gen double diff_hi_s = diff_hi     * `scale'

        quietly summarize diff_lo_s
        local ymin = r(min)
        quietly summarize diff_hi_s
        local ymax = r(max)
        local maxabs = max(abs(`ymin'), abs(`ymax'))
        if `maxabs'==. local maxabs = 0.1
        local pad = `maxabs'*0.20
        local ymin2 = -(`maxabs' + `pad')
        local ymax2 =  (`maxabs' + `pad')

        twoway ///
            (rarea diff_lo_s diff_hi_s p, sort fcolor(ebblue%20) lcolor(ebblue%50)) ///
            (line  diff_s    p, lcolor(navy)) ///
            , ///
              xscale(range(5 95)) ///
              yscale(range(`ymin2' `ymax2')) ///
              yline(0, lpattern(shortdash) lcolor(gs10)) ///
              ytitle("`=ustrunescape("\u0394(\u03B8) = W\u0302(d_\u03B8) \u2212 Y\u0304_RCT")'", size(medium)) ///
              xtitle("Percentile threshold of score") ///
              title("Policy value vs RCT — `outcomevar'", size(medlarge)) ///
              legend(order(2 1) rows(1) keygap(*0.8) ///
                     label(2 "TR − RCT") ///
                     label(1 "95% bootstrap band")) ///
              note("Shaded band: 95% stratified bootstrap percentile CI.", size(small)) ///
              graphregion(color(white)) ///
              name(G_bot_`outcomevar', replace) nodraw

        graph combine G_top_`outcomevar' G_bot_`outcomevar', ///
            rows(2) cols(1) xcommon imargin(small) ///
            ysize(5in) xsize(4in) ///
            name(G_`outcomevar', replace)

        graph export "`OUTDIR'/combined_`outcomevar'.png", replace width(`EX_W') height(`EX_H')
        graph export "`OUTDIR'/combined_`outcomevar'.pdf", replace
        graph save   "`OUTDIR'/combined_`outcomevar'.gph", replace
    }

    * Clean temps
    cap erase "`tmpplot'"
    cap erase "`bootrows'"
}

cap log close _all
