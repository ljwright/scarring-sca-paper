// Cell 1
clear
set more off
cd "C:/Users/liamj/Documents/Datasets/Next Steps - Wave 1- 8 - 2021-01-14"

global stata_fld  stata/stata13/eul
global dta_fld    Projects/Robustness of Scarring Effects/Submissions/Data
global act_fld    Projects/Activity Histories/Data

capture program drop prog_droplbl
program define prog_droplbl
    syntax anything
    
    qui foreach label of local anything {
        label list `label'
        if `r(min)'<0 {
            forval i=`r(min)'/-1 {
                local lbl: label `label' `i'
                if "`lbl'"!="`i'" label define `label' `i' "", modify 
            }
        }
    }
end

capture program drop prog_scale
program define prog_scale
	syntax varlist
	
	qui foreach var of local varlist{
		sum `var'
		replace `var' = (`var' - r(mean))/r(sd)
	}	
end

capture program drop prog_total
program define prog_total
	args newvar vlist
	
	egen `newvar' = rowtotal(`vlist')
	egen missing = rowmiss(`vlist')
	replace `newvar' = . if missing > 0
	drop missing
end

// Cell 2
use NSID W8DHANVQH using "${stata_fld}/ns8_2015_derived.dta", clear 
merge 1:1 NSID using "${stata_fld}/ns8_2015_main_interview.dta", ///
	nogen keepusing(W8FINWT W8HEIGHT W8PATIENCE W8CMSEX W8DACTIVITY) 
merge 1:1 NSID using "${stata_fld}/ns8_2015_self_completion.dta", ///
	nogen keepusing(W8GHQ12_* W8LOCUS0?) 
merge 1:1 NSID using "${stata_fld}/wave_five_lsype_young_person_2020.dta", ///
	nogen keepusing(W5SexYP) 
merge 1:1 NSID using "${stata_fld}/wave_four_lsype_young_person_2020.dta", ///
	nogen keepusing(W4ActivYP W4Boost W4ConcenYP W4DecideYP W4DepressYP ///
	W4DifficYP W4HappyYP W4Hea1CYP W4NoConfYP W4NoSleepYP W4ProbsYP ///
	W4SexYP W4StrainYP W4UsefulYP W4WthlessYP W4ethgrpYP W4Weight_MAIN_BOOST) 
merge 1:1 NSID using "${stata_fld}/wave_three_lsype_young_person_2020.dta", ///
	nogen keepusing(W3*bulrc W3pbulrc W3risk W3sexYP W3yschat1 ///
	W3yys1YP-W3yys12YP) 
merge 1:1 NSID using "${stata_fld}/wave_three_lsype_family_background_2020", ///
	nogen keepusing(W3cnssecfam W3hous12HH) 
merge 1:1 NSID using "${stata_fld}/wave_two_lsype_family_background_2020.dta", ///
	nogen keepusing(IMDRSCORE W2nssecfam W2Hous12HH) 
merge 1:1 NSID using "${stata_fld}/wave_two_lsype_young_person_2020.dta", ///
	nogen keepusing(W2Fat1YP W2Fat2YP W2Fat4YP W2Fat5YP W2Fat7YP W2Fat8YP ///
	W2FinWt W2SexYP W2YYS1YP-W2YYS12YP W2bulrc W2disabYP W2ethgrpYP ///
	W2ghq12scr W2hea1cYP W2pbulrc W2risk W2yschat1 W2concenYP-W2happyYP) 
merge 1:1 NSID using "${stata_fld}/wave_one_lsype_family_background_2020.dta", ///
	nogen keepusing(SampPSU SampStratum W1FinWt W1hiqualgMP W1hiqualgSP W1nssecfam W1hous12HH) 
merge 1:1 NSID using "${stata_fld}/wave_one_lsype_young_person_2020.dta", ///
	nogen keepusing(W1*bulrc W1disabYP W1ethgrpYP W1risk W1sexYP W1yschat1 ///
	W1yys1YP-W1yys12YP)
save "${dta_fld}/Collected Variables", replace


*  2. CLEAN DATA
/**/
use "${dta_fld}/Collected Variables", clear
numlabel, add

* Define Labels
label copy W1nssecfam NSSEC8
label define NSSEC8 -999 "" -99 "" -94 "" -92 "" -91 "", modify
label define Binary 0 "No" 1 "Yes"
label define Tenure 1 "Own House" 2 "Mortgage" 3 "Rent Council" 4 "Rent Private"
label define Disabled 1 "No" 2 "Yes, school not affected" 3 "Yes, school affected"
label define Ethnicity 1 "White" 2 "Mixed" 3 "Indian" 4 "Pakistani" 5 "Bangladeshi" /*
	*/ 6 "Black African" 7 "Black Caribbean" 8 "Other"
label define Education 1 "NVQ 5" 2 "NVQ 4" 3 "NVQ 3" 4 "NVQ 2" 5 "NVQ 1" 6 "No/Other Qual"
label define GenHealth 1 "Very Good" 2 "Fairly Good" 3 "Not Very Good" 4 "Not Good at All"
label define NSSEC3 1 "Higher" 2 "Intermediate" 3 "Routine" 4 "LTU"
label define LOC_Cat 1 "Neither" 2 "External" 3 "Internal"
label define Quintiles 1 "20th" 2 "40th" 3 "60th" 4 "80th" 5 "100th"
label define Bullied_Wave 0 "0" 1 "1" 2 "2" 3 "3"
label define Status_W8 1 "Employed" 2 "Education" 3 "Inactive" 4 "Unemployed"
label define Female 0 "Male" 1 "Female"

* Waves 1-8: Fixed Characteristics
egen Female=rowmax(W2SexYP W3sexYP W1sexYP W4SexYP W5SexYP)
replace Female=W8CMSEX if !inlist(Female,1,2)
replace Female=cond(Female>0,Female - 1,.)
egen Ethnicity=rowmax(W2ethgrpYP W1ethgrpYP W4ethgrpYP)
replace Ethnicity=cond(Ethnicity>0,Ethnicity,1)

label values Female Female
label values Ethnicity Ethnicity
drop *sexYP *SexYP *ethgrpYP W8CMSEX


* Wave 1 & 2: Common Time-Varying Variables
tab1 *disabYP
forval i=1/2{
	gen Disabled_W`i'=4-W`i'disabYP if inrange(W`i'disabYP,1,3)
	}
gen Disabled=max(Disabled_W1, Disabled_W2)
drop Disabled_W?
label values Disabled* Disabled
drop W*disabYP 

* Wave 2, 4 & 8: GHQ-12 & General Health
gen Survey_Weight_W2 = W2FinWt  if W2FinWt>0
gen Survey_Weight_W4 = W4Weight_MAIN_BOOST  if W4Weight_MAIN_BOOST>0
gen Survey_Weight_W8 = W8FINWT  if W8FINWT>0
drop W2FinWt W8FINWT W4Weight_MAIN_BOOST

gen GenHealth_W2=W2hea1cYP-2 if inrange(W2hea1cYP,3,6)
gen GenHealth_W4=W4Hea1CYP if inrange(W4Hea1CYP,1,4)
label values GenHealth* GenHealth
drop W2hea1cYP W4Hea1CYP

local ghq_positive 1 3 4 7 8 12
local ghq_negative 2 5 6 9 10 11

local w2ghq		W2concenYP W2nosleepYP W2usefulYP W2decideYP ///
				W2strainYP W2difficYP W2activYP W2probsYP ///
				W2depressYP W2noconfYP W2wthlessYP W2happyYP
local w4ghq 	W4ConcenYP W4NoSleepYP W4UsefulYP W4DecideYP ///
				W4StrainYP W4DifficYP W4ActivYP W4ProbsYP ///
				W4DepressYP W4NoConfYP W4WthlessYP W4HappyYP
local w8ghq		W8GHQ12_1 W8GHQ12_2 W8GHQ12_3 W8GHQ12_4 ///
				W8GHQ12_5 W8GHQ12_6 W8GHQ12_7 W8GHQ12_8 ///
				W8GHQ12_9 W8GHQ12_10 W8GHQ12_11 W8GHQ12_12

foreach w in 2 4 8{
	local i = 0
	foreach var of local w`w'ghq{
		local i = `i' + 1
		gen W`w'_GHQ_Likert_`i' = `var' - 1 if inrange(`var', 1, 4)
		gen W`w'_GHQ_Caseness_`i' = inrange(`var', 3, 4) if inrange(`var', 1, 4)
		if `w'==2	replace W`w'_GHQ_Caseness_`i' = 0 if `var'==-1
		
		gen W`w'_GHQ_Corrected_`i' = W`w'_GHQ_Caseness_`i'
	}
	
	foreach n of local ghq_negative{
		replace W`w'_GHQ_Corrected_`n' = 1 if W`w'_GHQ_Likert_`n'==1
	}
	
	foreach type in Likert Caseness Corrected{
		prog_total GHQ_W`w'_`type' W`w'_GHQ_`type'*
	}
}	

drop GHQ_W2*Likert W2ghq12scr `w2ghq' `w4ghq' `w8ghq'

* Waves 1, 2 and 3: Risk, Bullying, School Attitude, NS-SEC, Tenure
rename W3cnssecfam W3nssecfam
rename W2Hous12HH W2hous12HH
rename W*yys*, lower
rename W*YYS*, lower
recode w*yys* (min/-1 = .)
forval i=1/3{
	gen SchoolAtt_W`i'=W`i'yschat1 if inrange(W`i'yschat1,0,48)		
	gen Risk_W`i'=W`i'risk if inrange(W`i'risk,0,8)
	
	gen Bullied_W`i'=(W`i'bulrc==1 | W`i'pbulrc==1) if !missing(W`i'bulrc,W`i'pbulrc)
	replace Bullied_W`i' = . if  W`i'bulrc<0 & W`i'pbulrc<0
	
	gen NSSEC8_W`i'=W`i'nssecfam if inrange(W`i'nssecfam,1,8)
	gen NSSEC3_W`i'=1 if inrange(NSSEC8_W`i',1,2)
	replace NSSEC3_W`i'=2 if inrange(NSSEC8_W`i',3,4)
	replace NSSEC3_W`i'=3 if inrange(NSSEC8_W`i',5,7)
	replace NSSEC3_W`i'=4 if inlist(NSSEC8_W`i',8)
	
	gen Tenure_W`i' = 1 if W`i'hous12HH == 1
	replace Tenure_W`i' = 2 if inrange(W`i'hous12HH, 2, 3)
	replace Tenure_W`i' = 3 if inrange(W`i'hous12HH, 4, 5)
	replace Tenure_W`i' = 4 if inrange(W`i'hous12HH, 6, 8)
	}
gen Bullied_Waves=Bullied_W1+Bullied_W2+Bullied_W3

factor SchoolAtt_W?
predict SchoolAtt_Factor
factor Risk_W*
predict Risk_Factor

label values Bullied_Waves Bullied_Waves
label values NSSEC3_W? NSSEC3
label values Tenure_W? Tenure

drop W*risk W*yschat1 w*yys* W*bulrc W*nssecfam NSSEC8_W? W*hous12HH Bullied_W?


* Wave 1
ds W1*

gen ParentEduc5_W1 = W1hiqualgMP ///
    if inrange(W1hiqualgMP, 1, 7) & W1hiqualgSP == -98
replace ParentEduc5_W1 = min(W1hiqualgMP, W1hiqualgSP) ///
    if inrange(W1hiqualgMP, 1, 7) & inrange(W1hiqualgSP, 1, 7)
replace ParentEduc5_W1 = 5 if inrange(ParentEduc5_W1, 5, 7)
label define ParentEduc5_W1 1 "Degree" 2 "Other HE" 3 "A-Level" 4 "GCSE A-C" 5 "Other/None"
label values ParentEduc5_W1 ParentEduc5_W1

rename W1FinWt Survey_Weight_W1
rename SampPSU Survey_PSU_W1
rename SampStratum Survey_Stratum_W1

mca NSSEC3_W1 Tenure_W1 ParentEduc5_W1
predict SES_W1

drop W1* Tenure_W2 Tenure_W3 NSSEC3_W2 NSSEC3_W3


* Wave 2
gen IMD_W2=IMDRSCORE if IMDRSCORE>=0
xtile IMD_W2_Quintile=IMD_W2 [pweight=Survey_Weight_W2], n(5)
label values IMD_W2_Quintile Quintiles
drop IMDRSCORE

local Int_LOC W2Fat1YP W2Fat5YP W2Fat8YP
local Ext_LOC W2Fat2YP W2Fat4YP W2Fat7YP
foreach type in Int Ext{
	local i = 0
	foreach var of local `type'_LOC{
		local i = `i' + 1
		if "`type'"=="Int"{
			gen Int_LOC_W2_Item`i' = 5 - `var' if inrange(`var', 1, 2)
			replace Int_LOC_W2_Item`i' = 2 if `var'==-1
			replace Int_LOC_W2_Item`i' = 4 - `var' if inrange(`var', 3, 4)
		}	
		else{
			gen Ext_LOC_W2_Item`i' = `var' if inrange(`var', 3, 4)
			replace Ext_LOC_W2_Item`i' = 2 if `var'==-1
			replace Ext_LOC_W2_Item`i' = `var' - 1 if inrange(`var', 1, 2)
		} 		
	}
}
drop `Int_LOC' `Ext_LOC'
label variable Int_LOC_W2_Item1 "Q1. If someone is not a success in life, it is usually their own fault"
label variable Int_LOC_W2_Item2 "Q4. I can pretty much decide what will happen in my life"
label variable Int_LOC_W2_Item3 "Q6. If you work hard at something, you'll usually succeed"
label variable Ext_LOC_W2_Item1 "Q2. Even if I do well in school, I'll have a hard time"
label variable Ext_LOC_W2_Item2 "Q3. People like me don't have much of a chance in life"
label variable Ext_LOC_W2_Item3 "Q5. How well you get on in this world is mostly a matter of luck"

prog_total Int_LOC_W2_Sum Int_LOC_W2_Item*
prog_total Ext_LOC_W2_Sum Ext_LOC_W2_Item*
prog_total LOC_W2_Sum *LOC_W2_Item*


preserve
	polychoric *LOC*Item* [pweight = Survey_Weight_W2]
	matrix c = r(R)
	local n = r(N)
	factormat c, n(`n') factors(2)
	rotate, promax
	clear
	svmat2 e(r_L), names(col) rnames(item)
	gen eigen_1 = e(Ev)[1,1]
	gen eigen_2 = e(Ev)[1,2]
	save "${dta_fld}/loc_factors", replace
restore

* Wave 4
gen Survey_Boost=(W4Boost==1)
label values Survey_Boost Binary
drop W4Boost

* Wave 8
gen Interviewed_W8=(!missing(Survey_Weight_W8))
gen Has_GHQ_W8=(!missing(GHQ_W8_Likert))
gen Education_W8=6-W8DHANVQH if inrange(W8DHANVQH,1,5)
replace Education_W8=6 if inlist(W8DHANVQH,95,96)

label values Interviewed_W8 Has_GHQ_W8 Binary
label values Education* Education	
drop W8DHANVQH

sum W8HEIGHT if W8HEIGHT>0, d
gen Height_W8=W8HEIGHT*100 if inrange(W8HEIGHT,r(p1),r(p99))
gen Patience_W8=W8PATIENCE if inrange(W8PATIENCE,0,10)
drop W8HEIGHT W8PATIENCE

gen Status_W8 = 1 if inrange(W8DACTIVITY, 1, 4) | inlist(W8DACTIVITY, 8, 11, 12)
replace Status_W8 = 2 if inrange(W8DACTIVITY, 5, 6)
replace Status_W8 = 3 if inlist(W8DACTIVITY, 9, 10, 13, 14)
replace Status_W8 = 4 if W8DACTIVITY==5
label values Status_W8 Status_W8
drop W8DACTIVITY

recode W8LOCUS0? (min/-1 = .)
replace W8LOCUS0C = 5 - W8LOCUS0C
polychoric W8LOCUS0? [pweight = Survey_Weight_W8]
matrix c = r(R)
local n = r(N)
factormat c, n(`n') factors(2)
predict LOC_W8_Factor
drop W8LOCUS0?
	
* Format
drop Survey_PSU_W1 Survey_Stratum_W1 Survey_Weight_W1
ds *, alpha
order NSID `r(varlist)'
sort NSID
compress
save "${dta_fld}/Dataset", replace

// Cell 4
capture program drop prog_cleanact
program define prog_cleanact
	args lb ub
	
	replace End_MY = floor(End_MY)
	replace Start_MY = floor(Start_MY)
	replace Start_MY = `lb' if Start_MY<`lb'
	drop if Start_MY > End_MY

	by NSID (Spell), sort: replace Spell = _n
	expand 2 if Spell==1, gen(expand)
	replace Spell = 0 if expand==1
	replace Activity = .m if expand==1
	replace End_MY = Start_MY if expand==1
	replace Start_MY = `lb' if expand==1
	drop if End_MY == Start_MY
	drop expand
	by NSID (Spell), sort: replace Spell = _n

	by NSID (Spell), sort: gen gap = Start_MY[_n+1]-End_MY if _n<_N
	by NSID (Spell), sort: replace gap = `ub'-End_MY ///
		if End_MY<`ub' & _n==_N
		
	expand 2 if gap > 0 & !missing(gap), gen(expand)
	by NSID (Spell expand), sort: replace Spell = _n
	replace Start_MY = End_MY if expand==1
	replace End_MY = Start_MY + gap if expand==1
	replace Activity = .m if expand==1
	drop expand gap

	by NSID (Spell), sort: gen new_status = Activity!=Activity[_n-1]
	by NSID (Spell), sort: gen status_spell = sum(new_status)
	by NSID status_spell (Spell), sort: replace End_MY = End_MY[_N]
	by NSID status_spell (Spell), sort: keep if _n==1
	by NSID (status_spell), sort: replace Spell = _n
	drop new_status status_spell	
	
end

use "${act_fld}/Activity Histories", clear
keep NSID Spell Activity Start_MY End_MY

prog_cleanact "ym(2006, 9)" "ym(2013, 9)"

qui forval lb_y = 2006/2012{
forval ub_y = 2007/2013{
if inrange(`ub_y'-`lb_y', 1, 4) {
	local dur = `ub_y' - `lb_y'
	
	gen lb = max(ym(`lb_y', 10), Start_MY)
	gen ub = min(ym(`ub_y', 10), End_MY)
	replace lb = . if lb>ub
	replace ub = . if missing(lb)
	format lb ub %tm

	gen duration = ub - lb if Activity==3
	by NSID (Spell), sort: egen Cont_duration = max(duration)
	by NSID (Spell), sort: egen Cumul_duration = sum(duration)
	by NSID (Spell), sort: gen missing = Activity==.m & !missing(lb, ub)
	by NSID (Spell), sort: egen has_missing = max(missing)

	foreach m in 3 6 9 12{
		if `m'<`dur'*12 foreach type in Cont Cumul{
			local var Unem_`m'm_Age_`type'_`dur'_`lb_y'_`ub_y'
			gen `var' = `type'_duration >= `m' & !missing(`type'_duration)
			replace `var' = . if `var'==0 & has_missing == 1
		}
	}
	
	drop lb ub duration Cont_duration Cumul_duration missing has_missing
}
}
}

keep NSID Unem*
duplicates drop  
merge 1:1 NSID using "${dta_fld}/Dataset", keep(match using) nogen
compress
save "${dta_fld}/Dataset", replace


// Cell 5
use NSID Spell Activity Start_MY End_MY using "${act_fld}/Activity Histories", clear
merge m:1 NSID using "${stata_fld}/ns8_2015_main_interview.dta", ///
	nogen keepusing(W8INTMTH W8INTYEAR) 
gen IntDate_MY = ym(W8INTYEAR, W8INTMTH)
format *MY %tm
drop W8*
drop if missing(IntDate_MY)

prog_cleanact "ym(2006, 9)" "IntDate_MY"

gen FTE = Start_MY if Spell == 1
replace FTE = End_MY if Spell == 1 & Activity==2
by NSID (Spell), sort: ///
	replace FTE = cond(Activity==2 & Start_MY<FTE[_n-1]+12, End_MY, FTE[_n-1]) ///
	if _n>1
by NSID (Spell), sort: egen FTE_MY = max(FTE)
format FTE_MY %tm
drop FTE

keep if Start_MY>=FTE_MY
gen missing = 1 if Activity==.m & End_MY<FTE_MY+12
by NSID (Spell), sort: egen has_missing = max(missing)
drop if has_missing == 1
drop missing has_missing
by NSID (Spell), sort: replace Spell = _n

qui forval dur = 1/4{
	gen ub = min(End_MY, FTE + `dur'*12)
	replace ub = . if Start_MY>ub
	format ub %tm

	gen duration = ub - Start_MY if Activity==3
	by NSID (Spell), sort: egen Cont_duration = max(duration)
	by NSID (Spell), sort: egen Cumul_duration = sum(duration)
	by NSID (Spell), sort: gen missing = Activity==.m & !missing(ub)
	by NSID (Spell), sort: egen has_missing = max(missing)

	foreach m in 3 6 9 12{
		if `m'<`dur'*12 foreach type in Cont Cumul{
			local var Unem_`m'm_FTE_`type'_`dur'
			gen `var' = `type'_duration >= `m' & !missing(`type'_duration)
			replace `var' = . if `var'==0 & has_missing == 1
			replace `var' = . if FTE + `dur'*12 > IntDate_MY - 24
		}
	}
	drop ub duration Cont_duration Cumul_duration missing has_missing
}

keep NSID Unem*
duplicates drop
merge 1:1 NSID using "${dta_fld}/Dataset", keep(match using) nogen
compress
save "${dta_fld}/Dataset", replace


// Cell 7
use NSID Spell Activity Start_MY End_MY using "${act_fld}/Activity Histories", clear
prog_cleanact "ym(2008, 10)" "ym(2010, 5)"

gen lb = max(ym(2008, 10), Start_MY)
gen ub = min(ym(2010, 5), End_MY)
replace lb = . if lb>ub
replace ub = . if missing(lb)
format lb ub %tm

gen duration = ub - lb if Activity==3
by NSID (Spell), sort: egen max_duration = max(duration)
by NSID (Spell), sort: gen missing = Activity==.m & !missing(lb, ub)
by NSID (Spell), sort: egen has_missing = max(missing)

gen Unem_6Months = max_duration >= 6 & !missing(max_duration)
replace Unem_6Months = . if Unem_6Months==0 & has_missing == 1
drop lb ub duration max_duration missing has_missing

keep NSID Unem*
duplicates drop
merge 1:1 NSID using "${dta_fld}/Dataset", keep(match using) nogen
compress
foreach m in 3 6 9 12{
	label define Unem_`m'Months 0 "<`m' Months Unemployment" 1 "`m'+ Months Unemployment"
	label values Unem_`m'* Unem_`m'Months 
}
save "${dta_fld}/Dataset", replace