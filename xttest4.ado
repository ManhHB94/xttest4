*	Version 1.2		15 Aug 2025 by Manh H. B. (hbmanh9492@gmail.com)
*	Kezdi (2003) test for heteroscedasticity in FE model
*	Version 1.1 Correct degrees of freedom of tests 
*				in the presence of perfect multicollinearity
*	Version 1.2 works without vech() matrix function

cap program drop xttest4
program define xttest4, rclass
	version 11.0
	
	syntax 
	
	if "`e(cmd)'" != "xtreg" {
        display in red as error "last estimates not xtreg"
        exit 301
    }
	
	if "`e(model)'" != "fe" {
		di in red "last estimates not xtreg, fe"
		exit 301
	}
	
	preserve
	
	*	Get xvar form cmdline	(or use -indeplist- command)
	local cmdline `e(cmdline)'
		
	* 		Step 1:  cmd và depvar
	gettoken cmd rest : cmdline     // remove "xtreg"
	gettoken depvar rest : rest     // remove depvar

	* 		Step 2:  if 
	local pos = strpos("`rest'", " if ")
	if (`pos' > 0) {
		local rest = substr("`rest'", 1, `pos' - 1)
	}
	else {
	*		Step 3: in
		local pos = strpos("`rest'", " in ")
		if (`pos' > 0) {
			local rest = substr("`rest'", 1, `pos' - 1)
		}
		else {
	*		Step 4: "["" (weight)
			local pos = strpos("`rest'", "[")
			if (`pos' > 0) {
				local rest = substr("`rest'", 1, `pos' - 1)
			}
			else {
	*		Step 5: ","  (option)
				local pos = strpos("`rest'", ",")
				if (`pos' > 0) {
					local rest = substr("`rest'", 1, `pos' - 1)
				}
			}
		}
	}
	
	fvrevar `rest'				// treat factor variables, time series operators
	local xvar "`r(varlist)'"
	
	*	Predict e_it
	tempvar e
	qui predict `e' if e(sample), e	
	
	*	Mark balanced subsample
	tempvar touse
	qui xtset
	if `e(Tcon)'==0 {
		qui xtbalance2 `r(panelvar)' `r(timevar)' `xvar' if e(sample) , ///
			gen(`touse') o(NT)
	}
	
	else {
		qui gen byte `touse'=e(sample)
	}
	
	qui keep if `touse'
	
	*	id, time, obs
	tempvar id time obs
	qui gen `id' = `e(ivar)' if `touse'
	sort `id' `r(timevar)'
	qui by `id': gen `time'=_n if `touse'
	sort `id' `time'
	qui gen `obs'=_n
	
	*	Get N, T of balanced subsample
	tempname NT N T K
	qui xtsum `e(depvar)' 
	scalar `NT'=r(N)
	scalar `N'=r(n)
	scalar `T'=r(Tbar)
	scalar `K'= e(rank)-1

	*	Predict e_it
	tempvar s2
	qui sum `e' , d
	scalar `s2'=r(Var)*(r(N)-1)/(`NT'-`N')

	*	Time-Demeaning
	sort `id' `time'
	local xvar_dm
	foreach var of varlist `xvar' {
		tempvar `var'_m `var'_dm 
		qui by `id': egen ``var'_m'=mean(`var') 
		qui gen ``var'_dm' = `var'-``var'_m' 
		local xvar_dm `xvar_dm' ``var'_dm'
	}
	
	capture matrix drop `Omega' `V0' `V1' `V2' `V3' `E'
	tempname Omega V0 V1 V2 V3 E
		
	*	Omega ========================
	capture drop _t*
	qui tab `time' if `touse', gen(_t)
	mat opaccum `Omega' = _t* , gr(`id') opvar(`e') nocons
	capture drop _t*
	
	mat `Omega' = `Omega'/`N'
	
		
	*	V0 ========================
	mat opaccum `V0' = `xvar_dm' , gr(`id') opvar(`e') nocons
	mat `V0' = `V0'/`N'
	
	*	V1 ========================
	mat glsa `V1' = `xvar_dm' , gr(`id') gl(`Omega') r(`time') nocons
	mat `V1' = `V1'/`N'

	*	V2 ========================
	sort `obs'
	mat opaccum `V2' = `xvar_dm' , gr(`obs') opvar(`e') nocons
	mat `V2'=`V2'*`T'/(`NT'-`N')
	
	*	V3 ========================
	qui mat ac `V3' = `xvar_dm' , abs(`id') nocons
	mat `V3' = `V3'*`s2'/`N'
	
	*	vj = vech(Vj) =============
	tempname v0 v1 v2 v3
	
	if c(version) >= 18 {
		forvalues i=0/3 {
			mat `v`i'' = vech(`V`i'')
		}	    
	}
	else {
		forvalues i=0/3 {
			vec_h `V`i''
			mat `v`i'' = r(v)
		}
	}	

	*	Cj ========================
	tempname m C1 C2 C3 X c1i c2i c3i M1 M2 M3

	scalar `m'=rowsof(`v1')	
	mat `C1' = J(`m',`m',0)
	mat `C2' = J(`m',`m',0)
	mat `C3' = J(`m',`m',0)

	if c(version) >= 18 {
		forvalues i=1/`=`N'' {
			mkmat `e' if `id'==`i', mat(`E')
			mkmat `xvar_dm' if `id'==`i', mat(`X')

			* C1
			mat `M1' = `X''*`E'*`E''*`X'-`X''*`Omega'*`X'
			mat `c1i' = vech(`M1')
			mat `C1' = `C1' + `c1i'*`c1i''
			
			* C2
			mat `M2' = `X''*`E'*`E''*`X'-`X''*diag(vecdiag(`E'*`E''))*`X'
			mat `c2i' = vech(`M2')
			mat `C2' = `C2' + `c2i'*`c2i''
				
			* C3
			mat `M3' = `X''*`E'*`E''*`X'-`X''*`s2'*`X'
			mat `c3i' = vech(`M3')
			mat `C3' = `C3' + `c3i'*`c3i''
				
		}	    
	}
	
	else {
		forvalues i=1/`=`N'' {
			mkmat `e' if `id'==`i', mat(`E')
			mkmat `xvar_dm' if `id'==`i', mat(`X')

			* C1
			mat `M1' = `X''*`E'*`E''*`X'-`X''*`Omega'*`X'
			vec_h `M1'
			mat `c1i' = r(v)
			mat `C1' = `C1' + `c1i'*`c1i''
			
			* C2
			mat `M2' = `X''*`E'*`E''*`X'-`X''*diag(vecdiag(`E'*`E''))*`X'
			vec_h `M2'
			mat `c2i' = r(v)
			mat `C2' = `C2' + `c2i'*`c2i''
				
			* C3
			mat `M3' = `X''*`E'*`E''*`X'-`X''*`s2'*`X'
			vec_h `M3'
			mat `c3i' = r(v)
			mat `C3' = `C3' + `c3i'*`c3i''
				
		}	    
	}
		
	mat `C1' = `C1'/`N'
	mat `C2' = `C2'/`N'
	mat `C3' = `C3'/`N'
	
	*	Test 
	tempname df h1 h2 h3 h1_p h2_p h3_p
	scalar `df'=`K'*(`K'+1)/2+1

	forvalues i=1/3 {
		tempname _h`i'
		mat `_h`i'' = `N'*(`v`i''-`v0')'*invsym(`C`i'')*(`v`i''-`v0')
		scalar `h`i'' = `_h`i''[1,1]
		scalar `h`i'_p' = 1-chi2(`df', `h`i'')
	}
	
			
	*	Display 
	di
	di as text "Kézdi test for heteroscedasticity in Fixed Effects Model"
	if `e(Tcon)'==0 {
		di in red "The test is performed on a balanced subsample with N = " ///
			`N' ", T = " `T'
	}
	di as text "{hline 60}"
	di as text " Test for   |     Statistic    df     P-value"
	di as text "{hline 60}"
	di as result " H1 vs. Ha    " 	_col(19) %9.3f `h1' 	///
								_col(30) %5.0f `df' 	///
								_col(37) %9.3f `h1_p'
	di as result " H2 vs. Ha    " 	_col(19) %9.3f `h2'   	///
								_col(30) %5.0f `df' 	///
								_col(37) %9.3f `h2_p'
	di as result " H3 vs. Ha    " 	_col(19) %9.3f `h3'   	///
								_col(30) %5.0f `df' 	///
								_col(37) %9.3f `h3_p'
	di as text "{hline 60}"
	di "H1: Cross-sectional homoskedasticity."
	di "H2: Serially uncorrelated: e_it, x_it or both."
	di "H3: Homoskedasticity and serially uncorrelated."
	di "Ha: Heteroskedasticity."
	di
	
	*	Return list
	ret scalar df1   = `df'
	ret scalar df2   = `df'
	ret scalar df3   = `df'
	ret scalar h1   = `h1'
	ret scalar h2   = `h2'
	ret scalar h3   = `h3'
	ret scalar h1_p   = `h1_p'
	ret scalar h2_p   = `h2_p'
	ret scalar h3_p   = `h3_p'
	
end


*	Define vec_h() for STATA versions not allow vech()

capture program drop vec_h
program define vec_h, rclass
	syntax name(name=matname)
	
	// Check if matrix exists
	capture matrix list `matname'
	if _rc {
	    di as error "Matrix `matname' not found"
		exit 198
	}

	// Get row/col number
	local n = rowsof(`matname')
	local k = colsof(`matname')

	// Check square matrix
	if `n' != `k' {
		di as error "Matrix must be square"
		exit 198
	}

	// Create empty column vector
	tempname v
	matrix `v' = J(`=(`n'*(`n'+1))/2', 1, .)

	// Get the lower triangle elements
	local idx = 1
	forvalues i = 1/`n' {
		forvalues j = 1/`i' {
			matrix `v'[`idx',1] = `matname'[`i',`j']
			local idx = `idx' + 1
		}
	}

    // Return in r()
    return matrix v = `v'
end

