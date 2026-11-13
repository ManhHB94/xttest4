/*============================================================================
	Test for heteroscedasticity in fixed-T panel data models
============================================================================*/

*	Version 1.6		12 Nov 2026 by Manh Hoang-Ba (hbmanh9492@gmail.com)

/* -----------------------------------------------------------------
	Kezdi (2003) LM test for fixed-T one-way fixed-effects models
------------------------------------------------------------------*/
*	Version 1.1 Correct degrees of freedom of tests 
*				in the presence of perfect multicollinearity
*	Version 1.2 works without vech() matrix function
*	Version 1.3 works with unbalanced data

/* -----------------------------------------------------------------
	Add regression-based tests for fixed-T one-way effects models
------------------------------------------------------------------*/
*	Version 1.4 precise Time-demeaning steps
*				use generalized inverse when covariance matrix is singular
*				not allow weight in xtreg, fe
*				add regression-based test for FEM and REM

*	Version 1.5 correctly compute `e2_i' in Regression-based test

/* -----------------------------------------------------------------
	Speed / memory release
------------------------------------------------------------------*/
*	Version 1.6 Speed and memory
*				- all heavy computations are done in Mata
*				- singularity of C_j is decided by rank() instead of det()==0
*				- Kezdi test: df = rank of C_j (= K(K+1)/2 when C_j
*				  is nonsingular), NOT K(K+1)/2+1 as printed in Kezdi (2004).
*                 Also correct when regressors are omitted, 
*				  N < K(K+1)/2, or C_j is null (H2 with T = 2)
*				- Kezdi test: C2 uses the same T_i/(T_i-1) scaling as V2, so the
*				  influence function of H2 has mean zero under H0.


cap program drop xttest4
program define xttest4, rclass
	version 11.0
	
	syntax [, kezdi]
	
	if "`e(cmd)'" != "xtreg" {
		di as err "last estimates not xtreg"
		exit 301
	}
	
	if "`e(model)'" != "fe" & "`e(model)'" != "re" {
		di as err "last estimates not xtreg, fe or xtreg, re"
		exit 301
	}
	
	if "`e(model)'"=="fe" {
		if "`e(wexp)'" != "" {
			di as err "weights not allowed"
			exit 101
		}
		local est = 1
	}
	else {
		local est = 2
		if "`kezdi'" != "" {
			di as err "option kezdi not allowed with xtreg, re"
			exit 198
		}
	} 

	*	Panel variable
	local ivar    "`e(ivar)'"
	local cmdline `"`e(cmdline)'"'
	local gvar "`ivar'"
	
	*	Time variable (only used to order observations inside panels)
	local tvar ""
	cap qui xtset
	if _rc == 0 {
		if "`r(panelvar)'" == "`ivar'" local tvar "`r(timevar)'"
	}

	tempvar touse
	qui gen byte `touse' = e(sample)

	*	Get xvar form cmdline
	gettoken cmd rest : cmdline     // remove "xtreg"
	gettoken depvar rest : rest     // remove depvar	
	_xttest4_getindvars `rest'		// keep indep. variables
	fvrevar `r(varlist)'				
	qui _rmcoll `r(varlist)' if `touse' , forcedrop
	local xvar "`r(varlist)'"

	*	Predict e_it (fixed-effects model only; the RE tests use OLS residuals)
	tempvar e
	if `est' == 1 {
		qui predict double `e' if `touse' , e
	}
	

/*==========================================================================
		1. Test for heteroscedasticity in fixed-effects models
==========================================================================*/	
	
	if `est'==1 & "`kezdi'" != "" {
	
/*==========================================================================
			1.1. Kezdi LM test (balanced and unbalanced panels)
==========================================================================*/	
		tempname R df1 df2 df3 h1 h2 h3 h1_p h2_p h3_p
		
		*	Mata returns: 1 balanced? | 2 N | 3 NT | 4-6 h1 h2 h3
		*				  7-9 status of C1 C2 C3 (0 nonsingular, 1 singular, 2 null)
		*				  10-12 rank of C1 C2 C3 (= degrees of freedom)
		mata: _xttest4_kezdi("`e'", "`xvar'", "`touse'", "`gvar'", "`tvar'", "`R'")
		
		local bal = `R'[1,1]
		
		forvalues i = 1/3 {
			local c = 3 + `i'
			scalar `h`i''   = `R'[1,`c']
			local c = 9 + `i'
			scalar `df`i''  = `R'[1,`c']
			scalar `h`i'_p' = chi2tail(`df`i'', `h`i'')
		}
		
		*	Format statistics
		forvalues i = 1/3 {
			if `i' > 1 | `bal' {
				if `h`i'' > 99999 {
					local h`i'_fmt `"%9.0g"'
				}
				else	local h`i'_fmt `"%9.3f"'	
			}
		}

		*	Display
		di
		di as text "Test for heteroscedasticity in fixed-T panel data models"	
		
		forvalues i = 1/3 {
			if `i' > 1 | `bal' {
				local c = 6 + `i'
				if `R'[1,`c'] == 1 {
					di as err "Covariance matrix C`i' is singular. Use generalized inverse."
				}
				if `R'[1,`c'] == 2 {
					di as err "H`i' is not defined: covariance matrix C`i' is null (e.g. H2 with T = 2)."
				}
			}
		}
		
		di as text "{hline 50}"
		di as text "Hypothesis         Statistic    df     P-value"
		di as text "{hline 50}"
		
		if `bal' {
			di as result " H1 vs. Ha    " 	_col(19) `h1_fmt' `h1' 	///
											_col(30) %5.0f `df1' 	///
											_col(37) %9.3f `h1_p'		
		}
		
		di as result " H2 vs. Ha    " 	_col(19) `h2_fmt' `h2'   	///
										_col(30) %5.0f `df2' 		///
										_col(37) %9.3f `h2_p'
		di as result " H3 vs. Ha    " 	_col(19) `h3_fmt' `h3'   	///
										_col(30) %5.0f `df3' 		///
										_col(37) %9.3f `h3_p'
		di as text "{hline 50}"
		if `bal' {
			di as text "H1: Cross-sectional homoskedasticity"
		}
		else {
			di as error "H1 is missing due to unbalanced data"
		}
		di as text "H2: Serially uncorrelated: e_it, x_it or both"
		di as text "H3: Homoskedasticity and serially uncorrelated"
		di as text "Ha: Heteroskedasticity"
		di
			
		*	Return list
		if `bal' {
			ret scalar df1   = `df1'
			ret scalar h1    = `h1'
			ret scalar h1_p  = `h1_p'
		}
		
		ret scalar df2   = `df2'
		ret scalar h2    = `h2'
		ret scalar h2_p  = `h2_p'
		
		ret scalar df3   = `df3'
		ret scalar h3    = `h3'
		ret scalar h3_p  = `h3_p'	
	
	}	// end kezdi test

	
/*==========================================================================
			1.2. FEM: Regression-based test (Juhl & Sosa-Escudero, 2014)
==========================================================================*/
	
	else if `est' == 1 {
		tempvar e2 e2_dm
		tempname hold lm1 lm2 df1 df2 lm1_p lm2_p
		
		qui gen  double `e2'    = `e'^2 if `touse'
		qui gen  double `e2_dm' = .
		
		local xvar_dm
		foreach var of local xvar {
			tempvar d
			qui gen double `d' = .
			local xvar_dm `xvar_dm' `d'
		}
		
		*	e2_i = mean of e2 within panel;  e2_dm = e2 - e2_i;  x_dm = x - x_i
		mata: _xttest4_demean("`e2' `xvar'", "`e2_dm' `xvar_dm'", "`gvar'", "`touse'")
		
		*	regress overwrites e(): hold the estimates and give them back afterwards
		_estimates hold `hold'
		capture noisily {
		
*			H1: Var(e_it|X, u_i) = sigma_it^2
			qui reg `e2' `xvar' if `touse'
			scalar `lm1' = e(N) * e(r2)
			scalar `df1' = e(df_m)

*			H2: Var(e_it|X, u_i) = sigma_i^2
			qui reg `e2_dm' `xvar_dm' if `touse'
			scalar `lm2' = e(N) * e(r2)
			scalar `df2' = e(df_m)
		}
		local rc = _rc
		_estimates unhold `hold'
		if `rc' exit `rc'
		
		scalar `lm1_p' = chi2tail(`df1', `lm1')
		scalar `lm2_p' = chi2tail(`df2', `lm2')
		
*		Show results

		forvalues i = 1/2 {
			if `lm`i'' > 99999 {
					local lm`i'_fmt `"%9.0g"'
			}
			else	local lm`i'_fmt `"%9.3f"'	
		}
		

		di
		di as text "Test for heteroscedasticity in fixed-T panel data models"
		di as text "Model:" _col(8) "Fixed-effects"
		di as text "H0:" 	_col(8) "Var(e_it | X_i, u_i) = " as result "sigma^2"
		di as text "{hline 50}"
		di as text _col(1) "H1" /*_col(15) "|"*/ 	_col(19) "Statistic" ///
								_col(33) "df" 	_col(39)  "P-value"
		di as text "{hline 50}"
		
		di as result "sigma_it^2    " 	_col(19) `lm1_fmt' `lm1'   	///
										_col(30) %5.0f `df1' 		///
										_col(37) %9.3f `lm1_p'
										
		di as result "sigma_i^2    " 	_col(19) `lm2_fmt' `lm2'   	///
										_col(30) %5.0f `df2' 		///
										_col(37) %9.3f `lm2_p'
		di as text "{hline 50}"	
		
*		return scalar
		ret scalar lm1_p	= `lm1_p'
		ret scalar df1   	= `df1'
		ret scalar lm1   	= `lm1'
		
		ret scalar lm2_p	= `lm2_p'
		ret scalar df2   	= `df2'
		ret scalar lm2   	= `lm2'
		
	}

	
/*==========================================================================
		2. REM: Regression-based test (Montes-Rojas & Sosa-Escudero, 2011)
==========================================================================*/
	
	else {
		tempvar w w_dm_sq w_i_dm_sq
		tempname hold lm1 lm2 lm3 df1 df2 df3 lm1_p lm2_p lm3_p
		
		qui gen double `w_dm_sq'   = .
		qui gen double `w_i_dm_sq' = .
		
		local xvar_i
		local xvar_dm
		foreach var of local xvar {
			tempvar d1 d2
			qui gen double `d1' = .
			qui gen double `d2' = .
			local xvar_i  `xvar_i'  `d1'
			local xvar_dm `xvar_dm' `d2'
		}
		
		*	regress overwrites e(): hold the estimates and give them back afterwards
		_estimates hold `hold'
		capture noisily {
		
			qui reg `depvar' `xvar'	if `touse'	// OLS
			qui predict double `w' if `touse' , r	
			
			*	w_dm_sq  = (w - w_i)^2
			*	w_i_dm_sq = w_i^2 - mean_i(w_dm_sq)/(T_i - 1)   (one row per panel)
			*	x_i, x_dm = panel means / deviations of x
			mata: _xttest4_re_prep("`w'", "`xvar'", "`xvar_dm'", "`xvar_i'", ///
					"`w_dm_sq'", "`w_i_dm_sq'", "`gvar'", "`touse'")
			
*			H01: Var(e_it|X) = sigma_e^2 | Var(u_i|X) = sigma_u^2
			qui reg `w_dm_sq' `xvar_dm' if `touse'
			scalar `lm1' = e(N) * e(r2)
			scalar `df1' = e(df_m)

*			H02: Var(u_i|X) = sigma_u^2 | Var(e_it|X) = sigma_e^2
			qui reg `w_i_dm_sq' `xvar_i' if `touse'
			scalar `lm2' = e(N) * e(r2)
			scalar `df2' = e(df_m)
		}
		local rc = _rc
		_estimates unhold `hold'
		if `rc' exit `rc'
		
		scalar `lm1_p' = chi2tail(`df1', `lm1')
		scalar `lm2_p' = chi2tail(`df2', `lm2')
		
*		H03: Var(u_i|X) = sigma_u^2, Var(e_it|X) = sigma_e^2	
		scalar `lm3'   = `lm1' + `lm2'
		scalar `df3'   = `df1' + `df2'
		scalar `lm3_p' = chi2tail(`df3', `lm3')
		
*		Show results
		forvalues i = 1/3 {
			if `lm`i'' > 99999 {
					local lm`i'_fmt `"%9.0g"'
			}
			else	local lm`i'_fmt `"%9.3f"'	
		}
		

		di
		di as text "Test for heteroscedasticity in fixed-T panel data models"
		di as text "Model:" _col(8) "Random-effects"
		di as text "H1:" 	_col(8) "Var(u_i + e_it | X_i) = " as result "sigma_it^2"
		di as text "{hline 53}"
		di as text "Null hypothesis (H0)" _col(25) "Statistic" ///
								_col(39) "df" 	_col(45)  "P-value"
		di as text "{hline 53}"
		
		di as result "sigma_e^2 | sigma_u^2" 	_col(25) `lm1_fmt' `lm1'   	///
										_col(36) %5.0f `df1' 		///
										_col(43) %9.3f `lm1_p'
										
		di as result "sigma_u^2 | sigma_e^2" 	_col(25) `lm2_fmt' `lm2'   	///
										_col(36) %5.0f `df2' 		///
										_col(43) %9.3f `lm2_p'
										
		di as result "sigma_e^2 & sigma_u^2" _col(25) `lm3_fmt' `lm3'   	///
										_col(36) %5.0f `df3' 		///
										_col(43) %9.3f `lm3_p'										
		di as text "{hline 53}"	
		
*		return scalar
		ret scalar lm1_p	= `lm1_p'
		ret scalar df1   	= `df1'
		ret scalar lm1   	= `lm1'
		
		ret scalar lm2_p	= `lm2_p'
		ret scalar df2   	= `df2'
		ret scalar lm2   	= `lm2'
		
		ret scalar lm3_p	= `lm3_p'
		ret scalar df3   	= `df3'
		ret scalar lm3   	= `lm3'			
		
	}
	
end


/*=====================================================================
	Define _xttest4_getindvars to get independent varlist from cmdline
=====================================================================*/

cap pro drop _xttest4_getindvars
pro def _xttest4_getindvars, rclass
	cap syntax varlist(ts fv) [if] [in] [aw pw fw] [, *]
	return local varlist "`varlist'"
end


/*=====================================================================
	Mata engine
=====================================================================*/

mata: mata set matastrict on

mata:

// ---------------------------------------------------------------------
// Panel structure of the estimation sample.
//   p   : permutation that orders the sample by (panel, time, obs)
//   info: panelsetup() of the ordered sample
//   pid : panel number (1..N) of each observation, in original order
// ---------------------------------------------------------------------
void xttest4_groups(string scalar gvar, string scalar tvar, string scalar touse,
				 real colvector p, real matrix info, real colvector pid)
{
	real matrix    id, tm
	real colvector ids, st
	real scalar    n
	
	st_view(id, ., gvar, touse)
	n = rows(id)
	if (tvar != "") {
		st_view(tm, ., tvar, touse)
		p = order((id, tm, (1::n)), (1, 2, 3))
	}
	else {
		p = order((id, (1::n)), (1, 2))
	}
	ids  = id[p, 1]
	info = panelsetup(ids, 1)
	st   = J(n, 1, 0)
	st[info[., 1]] = J(rows(info), 1, 1)
	pid  = J(n, 1, .)
	pid[p] = runningsum(st)
}

// panel means of a vector given in original order -> N x 1
real colvector xttest4_pmean(real colvector v, real colvector p, real matrix info)
{
	return(panelsum(v[p], info) :/ (info[., 2] - info[., 1] :+ 1))
}

// vech(): lower triangle stacked by columns (same order as the loops below)
real colvector xttest4_vech(real matrix M)
{
	real scalar    j, k, c, d
	real colvector v
	
	d = rows(M)
	v = J(d*(d+1)/2, 1, .)
	c = 0
	for (j = 1; j <= d; j++) {
		for (k = j; k <= d; k++) {
			c = c + 1
			v[c] = M[k, j]
		}
	}
	return(v)
}

// in-place within-panel demeaning of X (rows already ordered by panel)
void xttest4_demean_sorted(real matrix X, real matrix info)
{
	real scalar    n, N, r1, r2, bs
	real colvector st, pid
	real matrix    M
	
	n  = rows(X)
	N  = rows(info)
	M  = panelsum(X, info) :/ (info[., 2] - info[., 1] :+ 1)
	st = J(n, 1, 0)
	st[info[., 1]] = J(N, 1, 1)
	pid = runningsum(st)
	bs = 100000
	for (r1 = 1; r1 <= n; r1 = r1 + bs) {
		r2 = min((r1 + bs - 1, n))
		X[|r1, 1 \ r2, .|] = X[|r1, 1 \ r2, .|] - M[pid[|r1 \ r2|], .]
	}
}

// Wald statistic N * d' C^{-1} d and its degrees of freedom = rank(C).
// Returns (h, status, rank):  status 0 = C nonsingular
//                             status 1 = C singular -> Moore-Penrose inverse
//                             status 2 = C numerically null -> h missing, rank 0
// rank = number of singular values above max(rows,cols)*epsilon(largest).
// -scale- is the size of the (unrestricted) second moment matrix of the
// influence functions; it is used to recognise a null C (e.g. H2 with T = 2,
// where V2 == V0 identically and C holds rounding noise only).
real rowvector xttest4_wald(real colvector d, real matrix C, real scalar N, real scalar scale)
{
	real colvector sv
	real scalar    tol, rk, h
	
	sv = svdsv(C)
	if (scale <= 0 | max(sv) <= 1e-20 * scale) return((., 2, 0))
	tol = max((rows(C), cols(C))) * epsilon(max(sv))
	rk  = sum(sv :> tol)
	if (rk < rows(C)) h = N * cross(d, pinv(C) * d)
	else              h = N * cross(d, invsym(C) * d)
	return((h, (rk < rows(C)), rk))
}

// ---------------------------------------------------------------------
// Kezdi (2003) LM test: balanced and unbalanced panels
// ---------------------------------------------------------------------
void _xttest4_kezdi(string scalar evar, string scalar xvars, string scalar touse,
					string scalar gvar, string scalar tvar, string scalar resname)
{
	real matrix    X, ev, Xs, info, Xb, ib, A, B1, B2, B3, Z, S, C1, C2, C3
	real matrix    Wb, Xk, Omega, Em, xe, CA
	real colvector p, pid, es, Ti, Tb, eb, e2b, wt, pr, v0, v1, v2, v3
	real rowvector w1, w2, w3, res
	real scalar    n, N, np, q, bal, T, s2, bs, i1, i2, a, b, nb, j, k, c, scale
	
	st_view(X,  ., tokens(xvars), touse)
	st_view(ev, ., evar, touse)
	n  = rows(X)
	np = cols(X)
	
	xttest4_groups(gvar, tvar, touse, p, info, pid)
	pid = J(0, 1, .)
	N   = rows(info)
	Ti  = info[., 2] - info[., 1] :+ 1
	
	//	ordered copies; X becomes time-demeaned in place (the only NT x K copy)
	Xs = X[p, .]
	es = ev[p, 1]
	p  = J(0, 1, .)
	xttest4_demean_sorted(Xs, info)
	
	bal = (min(Ti) == max(Ti))
	s2  = quadsum((es :- mean(es)):^2) / (n - N)
	
	if (bal) {
		T     = Ti[1]
		Em    = colshape(es, T)				// N x T,  row i = e_i
		Omega = quadcross(Em, Em) / N
	}
	
	q  = np*(np+1)/2
	v0 = J(q, 1, 0)
	v1 = J(q, 1, 0)
	v2 = J(q, 1, 0)
	C1 = J(q, q, 0)
	C2 = J(q, q, 0)
	C3 = J(q, q, 0)
	CA = J(q, q, 0)
	
	//	panels are processed in blocks so that memory stays bounded
	bs = max((1, floor(1000000/q)))
	for (i1 = 1; i1 <= N; i1 = i1 + bs) {
		i2  = min((i1 + bs - 1, N))
		nb  = i2 - i1 + 1
		a   = info[i1, 1]
		b   = info[i2, 2]
		ib  = info[|i1, 1 \ i2, 2|] :- (a - 1)
		Xb  = Xs[|a, 1 \ b, .|]
		eb  = es[|a \ b|]
		e2b = eb:^2
		Tb  = Ti[|i1 \ i2|]
		xe  = panelsum(Xb :* eb, ib)			// X_i'e_i
		
		A  = J(nb, q, .)
		B2 = J(nb, q, .)
		B3 = J(nb, q, .)
		if (bal) {
			B1 = J(nb, q, .)
			Wb = J(rows(Xb), np, .)
			for (k = 1; k <= np; k++) {
				Xk = colshape(Xb[., k], T)
				Wb[., k] = colshape(Xk * Omega, 1)		// Omega * x_ik
			}
		}
		
		//	column c of A, B1, B2, B3 holds element (k,j) of the K x K matrices
		//	  A  = X_i'e_i e_i'X_i
		//	  B1 = X_i' Omega X_i
		//	  B2 = X_i' diag(e_i^2) X_i
		//	  B3 = s2 * X_i'X_i
		c = 0
		for (j = 1; j <= np; j++) {
			for (k = j; k <= np; k++) {
				c = c + 1
				A[., c] = xe[., j] :* xe[., k]
				pr = Xb[., j] :* Xb[., k]
				S  = panelsum((pr, pr :* e2b), ib)
				B3[., c] = s2 :* S[., 1]
				B2[., c] = S[., 2]
				if (bal) B1[., c] = panelsum(Xb[., j] :* Wb[., k], ib)
			}
		}
		
		wt = editmissing(Tb :/ (Tb :- 1), 0)		// T_i/(T_i-1); 0 for singletons
		v0 = v0 + colsum(A)'
		v2 = v2 + colsum(B2 :* wt)'
		CA = CA + cross(A, A)
		
		//	influence function of H2: x'ee'x - T_i/(T_i-1) * sum_t e_t^2 x_t x_t'
		Z  = A - B2 :* wt
		C2 = C2 + cross(Z, Z)
		Z  = A - B3
		C3 = C3 + cross(Z, Z)
		if (bal) {
			v1 = v1 + colsum(B1)'
			Z  = A - B1
			C1 = C1 + cross(Z, Z)
		}
	}
	
	v0 = v0 / N
	v1 = v1 / N
	v2 = v2 / N
	v3 = xttest4_vech(quadcross(Xs, Xs) * s2 / N)
	C1 = C1 / N
	C2 = C2 / N
	C3 = C3 / N
	
	scale = max(diagonal(CA)) / N
	
	w1 = (., ., .)
	if (bal) w1 = xttest4_wald(v1 - v0, C1, N, scale)
	w2 = xttest4_wald(v2 - v0, C2, N, scale)
	w3 = xttest4_wald(v3 - v0, C3, N, scale)
	
	res = (bal, N, n, w1[1], w2[1], w3[1], w1[2], w2[2], w3[2], w1[3], w2[3], w3[3])
	st_matrix(resname, res)
}

// ---------------------------------------------------------------------
// Within-panel demeaning of variables:  dst_j = src_j - panel mean(src_j)
// ---------------------------------------------------------------------
void _xttest4_demean(string scalar src, string scalar dst,
					 string scalar gvar, string scalar touse)
{
	real colvector p, pid, v, mn
	real matrix    info
	string rowvector s, d
	real scalar    j
	
	s = tokens(src)
	d = tokens(dst)
	xttest4_groups(gvar, "", touse, p, info, pid)
	for (j = 1; j <= cols(s); j++) {
		v  = st_data(., s[j], touse)
		mn = xttest4_pmean(v, p, info)
		st_store(., d[j], touse, v - mn[pid])
	}
}

// ---------------------------------------------------------------------
// Variables of the Montes-Rojas & Sosa-Escudero regression-based test
// ---------------------------------------------------------------------
void _xttest4_re_prep(string scalar wvar, string scalar xvars, string scalar xdm,
					  string scalar xinames, string scalar wdmsq, string scalar widmsq,
					  string scalar gvar, string scalar touse)
{
	real colvector p, pid, Ti, first, w, wi, wsq, wsqi, col, x, xm
	real matrix    info
	string rowvector xs, xd, xi
	real scalar    n, j
	
	xttest4_groups(gvar, "", touse, p, info, pid)
	n     = rows(pid)
	Ti    = info[., 2] - info[., 1] :+ 1
	first = p[info[., 1]]					// one observation per panel
	
	w   = st_data(., wvar, touse)
	wi  = xttest4_pmean(w, p, info)			// w_i
	wsq = (w - wi[pid]):^2					// w_dm^2
	st_store(., wdmsq, touse, wsq)
	
	wsqi = xttest4_pmean(wsq, p, info)			// mean_i(w_dm^2)
	col  = J(n, 1, .)
	col[first] = wi:^2 - wsqi :/ (Ti :- 1)	// missing when T_i == 1
	st_store(., widmsq, touse, col)
	
	xs = tokens(xvars)
	xd = tokens(xdm)
	xi = tokens(xinames)
	for (j = 1; j <= cols(xs); j++) {
		x  = st_data(., xs[j], touse)
		xm = xttest4_pmean(x, p, info)
		st_store(., xd[j], touse, x - xm[pid])
		col = J(n, 1, .)
		col[first] = xm
		st_store(., xi[j], touse, col)
	}
}

end

mata: mata set matastrict off
