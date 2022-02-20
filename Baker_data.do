use "https://github.com/scunning1975/mixtape/raw/master/baker.dta", clear 

rename id i
rename year t
rename y Y
rename treat_date Ei
gen treated = Ei!=0 // create a treatment variable
gen never_tr = Ei==2004

gen K = t-Ei
gen D = K>=0 & Ei!=.

gen gvar = cond(Ei==., 0, Ei) // csdid: replace Ei==. w/ Ei==0

global ep event_plot
global g0 "default_look"
// global g1 xla(-$pre (1) $post) /*global g1 xla(-5(1)5)*/
global g2 xt("Periods since the event")
global g3 yt("Average causal effect")
global g  $g2 $g3
global t "together"

//  AZ: dd25-dd43
tsset i t
xtline Y, overlay legend(off) name(gY, replace)

// TWFE OLS 
reghdfe Y dd1 - dd23 dd25-dd48 , a(i t) cluster(state) /* dd24 is the ref year */
estimates store ols // saving the estimates for later
event_plot ols, stub_lag(dd#) $t $g0 graph_opt($g ti("OLS") name(gOLS, replace))

// Estimation with eventstudyinteract of Sun and Abraham (2020)
eventstudyinteract Y dd1 - dd23 dd25-dd48, vce(cluster state) absorb(i t) cohort(Ei) control_cohort(never_tr)
$ep e(b_iw)#e(V_iw), stub_lag(dd#) $t $g0 graph_opt($g ti("SA 20")  name(gSA, replace)) 
matrix sa_b = e(b_iw) // storing the estimates for later
matrix sa_v = e(V_iw)

// Estimation with did_imputation of Borusyak et al. (2021)
did_imputation Y i t Ei, horizons(0/23) pretrend(24) minn(0) autosample /**/
estimates store bjs // storing the estimates for later
$ep bjs, $t $g0 graph_opt($g ti("BJS 21") name(gBJS, replace))

// Estimation with did_multiplegt of de Chaisemartin and D'Haultfoeuille (2020)
did_multiplegt Y i t D, robust_dynamic dynamic(23) placebo(24) breps(20) cluster(state) 
event_plot e(estimates)#e(variances), stub_lag(Effect_#) stub_lead(Placebo_#) $t $g0 graph_opt($g ti("CD 20") name(gCD, replace))
matrix dcdh_b = e(estimates) // storing the estimates for later
matrix dcdh_v = e(variances)

// Estimation with csdid of Callaway and Sant'Anna (2020)
csdid Y, ivar(i) time(t) gvar(gvar) notyet
estat event, estore(cs) // this produces and stores the estimates at the same time
$ep cs, stub_lag(Tp#) stub_lead(Tm#) $t $g0 graph_opt($g ti("CS 20") name(gCS, replace))

/* GB: bacondecomp */
set matsize 11000
bacondecomp Y D, ddetail legend(off) name(gGB, replace)

/* did2s (Gardner 2021) */	
did2s Y, first_stage(i t) second_stage(dd1 - dd23 dd25-dd48) treatment(D) cluster(state)
$ep, stub_lag(dd#) $t $g0 graph_opt($g ti("Gardner 21") name(gG, replace)) 
// did2s Y, first_stage(i t) second_stage(F_* L_*) treatment(D) cluster(i)
// $ep, stub_lag(L_#) stub_lead(F_#) $t $g0 graph_opt($g ti("Gardner 21") name(gG, replace)) 
matrix did2s_b = e(b)
matrix did2s_v = e(V)

/* stackedev (Cengiz, Dube, Lindner, Zipperer 2019) */
rename dd24 ref  // reference year
stackedev Y dd* ref, cohort(Ei) time(t) never_treat(never_tr) unit_fe(i) clust_unit(state)
$ep, stub_lag(dd#) $t $g0 graph_opt($g ti("CDLZ 19") name(gCDLZ, replace)) 
matrix stackedev_b = e(b)
matrix stackedev_v = e(V)
// stackedev Y F_* L_* ref, cohort(Ei) time(t) never_treat(never_tr) unit_fe(i) clust_unit(i)
// $ep, stub_lag(L_#) stub_lead(F_#) $t $g0 graph_opt($g ti("CDLZ 19") name(gCDLZ, replace)) 

/* gY gBJS gCD gCS gSA gOLS gGB gG gCDLZ */
graph combine gY gOLS gGB gBJS gCD gCS gSA gG gCDLZ, ycommon name(combined, replace)
