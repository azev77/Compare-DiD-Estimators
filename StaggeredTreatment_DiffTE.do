/*https://github.com/borusyak/did_imputation/blob/main/five_estimators_example.do*/
// Generate a complete panel of 300 units observed in 15 periods
clear all
timer clear
set seed 10
global T = 15
global I = 100   /*300*/

global pre  5 /*6*/
global post 8


global ep event_plot
global g0 "default_look"
global g1 xla(-$pre (1) $post) /*global g1 xla(-5(1)5)*/
global g2 xt("Periods since the event")
global g3 yt("Average causal effect")
global g  $g1 $g2 $g3
global t "together"

set obs `=$I*$T'
gen i = int((_n-1)/$T )+1 					// unit id
gen t = mod((_n-1),$T )+1					// calendar period
tsset i t

// Randomly generate treatment rollout years uniformly across Ei=10..16 (note that periods t>=16 would not be useful since all units are treated by then)
gen     Ei = 7 if t==1 & i>=1 & i <=35			// year when unit is first treated
replace Ei = 11 if t==1 & i>=36 & i <=70			// year when unit is first treated
bys i (t): replace Ei = Ei[1]
gen K = t-Ei 								// "relative time", i.e. the number periods since treated (could be missing if never-treated)
gen D = K>=0 & Ei!=. 						// treatment indicator

// Generate the outcome with parallel trends and heterogeneous treatment effects
gen tau = 0                                    // TE
replace tau = 20 if  D==1 & i>=1  & i <=35
replace tau = 60 if  D==1 & i>=36 & i <=70
gen eps = rnormal()							// error term
// gen Y = i + 3*t + tau*D + 3*eps 				// the outcome (FEs play no role since all methods control for them)
gen Y = i/100 + tau*D + 3*eps

// save five_estimators_data, replace

gen gvar = cond(Ei==., 0, Ei) // csdid: replace Ei==. w/ Ei==0
//  sum Ei
gen never_tr = Ei==. // never_tr =1 for never-treated cohort

/* Gen leads & lags of treatment */
forval x = 1/$pre {  // drop the first lead
	gen     F_`x' = K == -`x'
	replace F_`x' = 0 if never_tr==1
}
forval x = 0/$post {
	gen     L_`x' = K ==  `x'
	replace L_`x' = 0 if never_tr==1
}
rename F_1 ref  // reference year

/* Previous Leads/Lags, AZ: I think above is correct
cap drop F_* // leads
cap drop L_*      //lags
forval x = 1/14 {  
	gen F_`x' = K == -`x'
}
forval x = 0/5 {
	gen L_`x' = K ==  `x'
}
rename F_1 ref  // reference year
*/


xtline Y, overlay legend(off) name(gY, replace)

/* BJS & CD don't want too many pre-periods*/

// Estimation with did_imputation of Borusyak et al. (2021)
did_imputation Y i t Ei, horizons(0/$post) pretrend($pre) minn(0) /**/
estimates store bjs // storing the estimates for later
$ep bjs, $t $g0 graph_opt($g ti("BJS 21") name(gBJS, replace))

// Estimation with did_multiplegt of de Chaisemartin and D'Haultfoeuille (2020)
did_multiplegt Y i t D, robust_dynamic dynamic($post) placebo($pre) breps(20) cluster(i) 
event_plot e(estimates)#e(variances), stub_lag(Effect_#) stub_lead(Placebo_#) $t $g0 graph_opt($g ti("CD 20") name(gCD, replace))
matrix dcdh_b = e(estimates) // storing the estimates for later
matrix dcdh_v = e(variances)

// Estimation with csdid of Callaway and Sant'Anna (2020)
csdid Y, ivar(i) time(t) gvar(gvar) notyet
estat event, estore(cs) // this produces and stores the estimates at the same time
$ep cs, stub_lag(Tp#) stub_lead(Tm#) $t $g0 graph_opt($g ti("CS 20") name(gCS, replace))

// Estimation with eventstudyinteract of Sun and Abraham (2020)
eventstudyinteract Y L_* F_*, vce(cluster i) absorb(i t) cohort(Ei) control_cohort(never_tr)
$ep e(b_iw)#e(V_iw), stub_lag(L_#) stub_lead(F_#) $t $g0 graph_opt($g ti("SA 20")  name(gSA, replace)) 
matrix sa_b = e(b_iw) // storing the estimates for later
matrix sa_v = e(V_iw)

// TWFE OLS estimation (which is correct here because of treatment effect homogeneity). Some groups could be binned.
reghdfe Y L_* F_*, a(i t) cluster(i) /* cluster(i t) cluster(t) */
estimates store ols // saving the estimates for later
$ep ols, stub_lag(L_#) stub_lead(F_#) $t $g0 graph_opt($g ti("OLS") name(gOLS, replace))  


// More stuff
/* GB: bacondecomp */
bacondecomp Y D, ddetail legend(off) name(gGB, replace)

/* did2s (Gardner 2021) */	
did2s Y, first_stage(i t) second_stage(F_* L_*) treatment(D) cluster(i)
$ep, stub_lag(L_#) stub_lead(F_#) $t $g0 graph_opt($g ti("Gardner 21") name(gG, replace)) 
matrix did2s_b = e(b)
matrix did2s_v = e(V)


/* stackedev (Cengiz, Dube, Lindner, Zipperer 2019) */
stackedev Y F_* L_* ref, cohort(Ei) time(t) never_treat(never_tr) unit_fe(i) clust_unit(i)
$ep, stub_lag(L_#) stub_lead(F_#) $t $g0 graph_opt($g ti("CDLZ 19") name(gCDLZ, replace)) 
matrix stackedev_b = e(b)
matrix stackedev_v = e(V)

/* gY gBJS gCD gCS gSA gOLS gGB gG gCDLZ */
graph combine gY gOLS gGB gBJS gCD gCS gSA gG gCDLZ, ycommon name(combined, replace)

// Construct the vector of true average treatment effects by the number of periods since treatment
matrix btrue = J(1,9,.)
matrix colnames btrue = tau0 tau1 tau2 tau3 tau4 tau5 tau6 tau7 tau8
qui forvalues h = 0/8 {
	sum tau if K==`h'
	matrix btrue[1,`h'+1]=r(mean)
}

// Combine all plots using the stored estimates
event_plot /// 
		      btrue# bjs  dcdh_b#dcdh_v cs  sa_b#sa_v  did2s_b#did2s_v stackedev_b#stackedev_v ols, ///
	stub_lag( tau#   tau# Effect_#      Tp# L_#        L_#             L_#                     L_#) ///
	stub_lead(pre#   pre# Placebo_#     Tm# F_#        F_#             F_#                     F_#) ///
	plottype(scatter) ciplottype(rcap) ///
	together perturb(-0.325(0.1)0.325) trimlead(5) noautolegend ///
	graph_opt(  ///
	title("Event study estimators in a simulated panel", size(med)) ///
	xtitle("Periods since the event", size(small)) ///
	ytitle("Average causal effect", size(small)) xlabel(-$pre(1)$post)  ///
	legend(order(1 "Truth" 2 "BJS" 4 "dCdH" ///
				6 "CS" 8 "SA" 10 "G" 12 "CDLZ" 14 "TWFE") rows(2) position(6) region(style(none))) ///
	/// the following lines replace default_look with something more elaborate
		xline(-0.5, lcolor(gs8) lpattern(dash)) yline(0, lcolor(gs8)) graphregion(color(white)) bgcolor(white) ylabel(, angle(horizontal)) ///
	) 	///
	lag_opt1(msymbol(+) color(black)) lag_ci_opt1(color(black)) ///
	lag_opt2(msymbol(O) color(cranberry)) lag_ci_opt2(color(cranberry)) ///
	lag_opt3(msymbol(Dh) color(navy)) lag_ci_opt3(color(navy)) ///
	lag_opt4(msymbol(Th) color(forest_green)) lag_ci_opt4(color(forest_green)) ///
	lag_opt5(msymbol(Sh) color(dkorange)) lag_ci_opt5(color(dkorange)) ///
	lag_opt6(msymbol(Th) color(blue)) lag_ci_opt6(color(blue)) ///
	lag_opt7(msymbol(Dh) color(red)) lag_ci_opt7(color(red))  ///
	lag_opt8(msymbol(Oh) color(purple)) lag_ci_opt8(color(purple))



event_plot    btrue# bjs  dcdh_b#dcdh_v cs  sa_b#sa_v ols, ///
	stub_lag( tau#   tau# Effect_#      Tp# L_#       L_#) ///
	stub_lead(pre#   pre# Placebo_#     Tm# F_#       F_#) plottype(scatter) ciplottype(rcap) ///
	 noautolegend /// together perturb(-0.325(0.13)0.325) trimlead(5)
	graph_opt(title("Event study estimators in a simulated panel (300 units, 15 periods)", size(medlarge)) ///
		xtitle("Periods since the event") ytitle("Average causal effect") xlabel(-$pre(1)$post) /// ylabel(0(1)3) 
		legend(order(1 "True value" 2 "Borusyak et al." 4 "de Chaisemartin-D'Haultfoeuille" ///
				6 "Callaway-Sant'Anna" 8 "Sun-Abraham" 10 "OLS") rows(3) region(style(none))) ///
	/// the following lines replace default_look with something more elaborate
		xline(-0.5, lcolor(gs8) lpattern(dash)) yline(0, lcolor(gs8)) graphregion(color(white)) bgcolor(white) ylabel(, angle(horizontal)) ///
	) ///
	lag_opt1(msymbol(+) color(cranberry)) lag_ci_opt1(color(cranberry)) ///
	lag_opt2(msymbol(O) color(cranberry)) lag_ci_opt2(color(cranberry)) ///
	lag_opt3(msymbol(Dh) color(navy)) lag_ci_opt3(color(navy)) ///
	lag_opt4(msymbol(Th) color(forest_green)) lag_ci_opt4(color(forest_green)) ///
	lag_opt5(msymbol(Sh) color(dkorange)) lag_ci_opt5(color(dkorange)) ///
	lag_opt6(msymbol(Oh) color(purple)) lag_ci_opt6(color(purple)) 
graph export "five_estimators_example.png", replace
