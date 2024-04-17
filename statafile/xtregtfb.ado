////////////////////////////////////////////////////////////////////////////////
// STATA command for Chen and Vogelsang (2023)
// An adaption from xtregtwo.ado for Chiang, Hansen, Sasaki (2022)
////////////////////////////////////////////////////////////////////////////////
//
//

!* version 17.0  17Apr2024
program define xtregtfb, eclass
    version 17.0
 
    syntax varlist(numeric) [if] [in] [, NOConstant FE lag(real -1) SE(string) level(real 0.05) bm(real 1000) rep(real 2000) whichvar(real 1)]
    marksample touse
 
	qui xtset
	local panelid = r(panelvar)
	local timeid  = r(timevar)

    gettoken depvar indepvars : varlist
    _fv_check_depvar `depvar'
    fvexpand `indepvars' 
    local cnames `r(varlist)'

    tempname b V N T NT Mhat cvDKA cvBCCHS

	local const = 0
	if "`noconstant'" == "" {
	  local const = 1
	}
	
	local fixedeffect = 1
	if "`fe'" == "" {
	  local fixedeffect = 0
	}
	
	//0-chs, 1-bcchs, 2-dka, 3-bcchsfb, 4-dkafb
	local type = 1
	if "`se'" == "chs" {
		local type = 0
	}
	if "`se'" == "bcchs" {
		local type = 1
	}
	if "`se'" == "dka" {
		local type = 2
	}
	if "`se'" == "bcchsfb" {
		local type = 3
	}
	if "`se'" == "dkafb" {
		local type = 4
	}
	if "`se'" == "ci" {
		local type = -1
	}
	if "`se'" == "ct" {
		local type = -2
	}
	if "`se'" == "ehw" {
		local type = -3
	}
	if "`se'" == "dk" {
		local type = -4
	}
	
	mata: estimation("`depvar'","`cnames'","`panelid'","`timeid'",		 ///
					 "`touse'","`b'","`V'","`N'","`T'","`NT'","`Mhat'",     ///
					 `const',`fixedeffect', `lag', `type', `level',      ///
					 "`cvBCCHS'", "`cvDKA'", `bm', `rep', `whichvar')
	
	local cnames `cnames'
	if "`noconstant'" == "" & "`fe'" == ""{
		local cnames "`cnames' _cons"
	}
	matrix colnames `b' = `cnames'
	matrix colnames `V' = `cnames'
	matrix rownames `V' = `cnames'
	
	ereturn post `b' `V', esample(`touse') buildfvinfo
	ereturn scalar N    = `N'
	ereturn scalar T    = `T'
	ereturn scalar NT   = `NT'
	ereturn scalar Mhat   = `Mhat'
	ereturn scalar cvBCCHS = `cvBCCHS'
	ereturn scalar cvDKA = `cvDKA'
	ereturn scalar lag = `Mhat'
	ereturn scalar level = `level'
	ereturn scalar bm = `bm'
	ereturn scalar rep = `rep'
	ereturn local  cmd  "xtregtfb"
	ereturn display
	di "cvBCCHS=cvDKA=" `cvDKA'
	di "Mhat=" `Mhat'
	di "Note: This command provides bias-correction through fixed-b approximation for the two-way clustering robust standard error."
	di "Reference: K. Chen and T.J. Vogelsang (2023) Fixed-b Asymptotics for Panel Models with Two-Way Clustering. Working Paper."
end
////////////////////////////////////////////////////////////////////////////////
 
mata:
//////////////////////////////////////////////////////////////////////////////// 
// Estimation Function
void estimation(string scalar depvar, 	string scalar indepvars, 			 ///
                string scalar panelid, 	string scalar timeid,    			 ///
				string scalar touse, 	string scalar bname,	 			 ///
				string scalar Vname,   	string scalar nname,	 			 ///
				string scalar tname,	string scalar ntname,				 ///
				string scalar mname,	real scalar constant,				 ///
				real scalar fe     ,    real scalar la,                      ///
				real scalar ty     ,    real scalar le,                      ///
				string scalar cvBCCHSname, string scalar cvDKAname,           ///
				real scalar bm     ,    real scalar rep,                   ///
				real scalar whichv )
{
    Y = st_data(., depvar, touse)
    X = st_data(., indepvars, touse)
    t = st_data(., timeid, touse)
	i = st_data(., panelid, touse)
	uniq_i = uniqrows(i)
	uniq_t = uniqrows(t)
	N = rows(uniq_i)
	T = rows(uniq_t)
	
	if( constant > 0.5 & fe < 0.5){
		ones = J(rows(X),1,1)
		X = X,ones
	}
	
//////////////////////////////////////////////////////////////////////////////// 
// Within-Transformation If Fixed Effect
	if( fe > 0.5 ){
		for( uidx = 1 ; uidx <= N ; uidx++ ){
			idx = uniq_i[uidx,1]
			for( jdx = 1 ; jdx <= cols(X) ; jdx++ ){
				X[select((1..rows(X))',i:==idx),jdx] = X[select((1..rows(X))',i:==idx),jdx] :- mean(X[select((1..rows(X))',i:==idx),jdx])
			}	
			Y[select((1..rows(X))',i:==idx),  1] = Y[select((1..rows(X))',i:==idx),  1] :- mean(Y[select((1..rows(X))',i:==idx),  1])
		}
		uniq_i = uniqrows(i)
		uniq_t = uniqrows(t)
		N = rows(uniq_i)
		T = rows(uniq_t)
	}

//////////////////////////////////////////////////////////////////////////////// 
// OLS and Residual
	beta = luinv(X'*X)*X'*Y	
	Uhat = Y :- X*beta
	
//////////////////////////////////////////////////////////////////////////////// 
// Variance Estimation

	// Omega Hat 1st Term (Arellano term)
	Omega_hat_1 = J(cols(X),cols(X),0)
	for( uidx = 1 ; uidx <= N ; uidx++ ){
		idx = uniq_i[uidx,1]
		for( jdx = 1 ; jdx <= cols(X) ; jdx++ ){
			for( kdx = 1 ; kdx <= cols(X) ; kdx++ ){
				Omega_hat_1[jdx,kdx] = Omega_hat_1[jdx,kdx] + sum((X[select((1..rows(X))',i:==idx),jdx]:*Uhat[select((1..rows(X))',i:==idx),1])) * sum((X[select((1..rows(X))',i:==idx),kdx]:*Uhat[select((1..rows(X))',i:==idx),1]))	
			}
		}
	} 
	
	// Omega Hat 2nd Term (leading term of DK)
	Omega_hat_2 = J(cols(X),cols(X),0)
	for( utdx = 1 ; utdx <= T ; utdx++ ){
		tdx = uniq_t[utdx,1]
		for( jdx = 1 ; jdx <= cols(X) ; jdx++ ){
			for( kdx = 1 ; kdx <= cols(X) ; kdx++ ){
				Omega_hat_2[jdx,kdx] = Omega_hat_2[jdx,kdx] + sum((X[select((1..rows(X))',t:==tdx),jdx]:*Uhat[select((1..rows(X))',t:==tdx),1])) * sum((X[select((1..rows(X))',t:==tdx),kdx]:*Uhat[select((1..rows(X))',t:==tdx),1]))
	}}} 
	
	// Omega Hat 3rd Term (leading term of NW)
	Omega_hat_3 = J(cols(X),cols(X),0)
	for( jdx = 1 ; jdx <= cols(X) ; jdx++ ){
		for( kdx = 1 ; kdx <= jdx ; kdx++ ){
			Omega_hat_3[jdx,kdx] = (X[,jdx]:*Uhat[,1])' * (X[,kdx]:*Uhat[,1])
			Omega_hat_3[kdx,jdx] = Omega_hat_3[jdx,kdx]
	}} 
	
	// Select M Hat
	if (fe == 1 | constant == 0){
		rho_hat = J(cols(X), 1, 0) //no constant term
		for( kdx = 1 ; kdx <= cols(X) ; kdx++ ){	
		utdx = 2
		tdx = uniq_t[utdx,1]
		tdx_lag = uniq_t[utdx-1,1]
		S = sum(X[select((1..rows(X))',t:==tdx),kdx]:*Uhat[select((1..rows(X))',t:==tdx),1]) 
		S_lag = sum(X[select((1..rows(X))',t:==tdx_lag),kdx]:*Uhat[select((1..rows(X))',t:==tdx_lag),1])
		for( utdx = 3 ; utdx <= T ; utdx++ ){
			tdx = uniq_t[utdx,1]
			tdx_lag = uniq_t[utdx-1,1]
			S = S \ sum(X[select((1..rows(X))',t:==tdx),kdx]:*Uhat[select((1..rows(X))',t:==tdx),1])
			S_lag = S_lag \ sum(X[select((1..rows(X))',t:==tdx_lag),kdx]:*Uhat[select((1..rows(X))',t:==tdx_lag),1])
		}
		rho_hat[kdx] = ( luinv((J(rows(S_lag),1,1),S_lag)' * (J(rows(S_lag),1,1),S_lag)) * (J(rows(S_lag),1,1),S_lag)'*S )[2,1] //with constant
		//rho_hat[kdx] = ( luinv((S_lag)' * (S_lag)) * (S_lag)'*S )[1,1] //without constant
		}
	}
	
	if (fe == 0 & constant == 1){
		rho_hat = J(cols(X) - 1, 1, 0) //no constant term
		for( kdx = 1 ; kdx <= cols(X) - 1 ; kdx++ ){	
		utdx = 2
		tdx = uniq_t[utdx,1]
		tdx_lag = uniq_t[utdx-1,1]
		S = sum(X[select((1..rows(X))',t:==tdx),kdx]:*Uhat[select((1..rows(X))',t:==tdx),1]) //no sum over units?
		S_lag = sum(X[select((1..rows(X))',t:==tdx_lag),kdx]:*Uhat[select((1..rows(X))',t:==tdx_lag),1])
		for( utdx = 3 ; utdx <= T ; utdx++ ){
			tdx = uniq_t[utdx,1]
			tdx_lag = uniq_t[utdx-1,1]
			S = S \ sum(X[select((1..rows(X))',t:==tdx),kdx]:*Uhat[select((1..rows(X))',t:==tdx),1])
			S_lag = S_lag \ sum(X[select((1..rows(X))',t:==tdx_lag),kdx]:*Uhat[select((1..rows(X))',t:==tdx_lag),1])
		}
		//rho_hat[kdx] = ( luinv((J(rows(S_lag),1,1),S_lag)' * (J(rows(S_lag),1,1),S_lag)) * (J(rows(S_lag),1,1),S_lag)'*S )[2,1] //with constant
		rho_hat[kdx] = ( luinv((S_lag)' * (S_lag)) * (S_lag)'*S )[1,1] //without constant
		}
	}
	
	
	mhat = floor( 1.8171 * ( sum(rho_hat:^2:/(1:-rho_hat):^4) / sum((1:-rho_hat:^2):^2:/(1:-rho_hat):^4) )^(1/3) * T^(1/3) ) + 1
	
	m=mhat  //m takes mhat in default unless lag() is specified.
	if (la != -1) {
		m=la
	}
	
	// Omega Hat 4th Term
	Omega_hat_4 = J(cols(X),cols(X),0)
	for( j = 1 ; j <= m-1 ; j++ ){
		temp_Omega_hat_4 = J(cols(X),cols(X),0)
		w = 1 - (j/m)
		for( utdx = (j+1) ; utdx <= T ; utdx++ ){
			tdx = uniq_t[utdx,1]
			tdxminus = uniq_t[utdx-j,1]
			for( jdx = 1 ; jdx <= cols(X) ; jdx++ ){
				for( kdx = 1 ; kdx <= cols(X) ; kdx++ ){
					temp_Omega_hat_4[jdx,kdx] = temp_Omega_hat_4[jdx,kdx] + sum((X[select((1..rows(X))',t:==tdx),jdx]:*Uhat[select((1..rows(X))',t:==tdx),1])) * sum((X[select((1..rows(X))',t:==tdxminus),kdx]:*Uhat[select((1..rows(X))',t:==tdxminus),1]))
		}}} 
		temp_Omega_hat_4 = w * temp_Omega_hat_4
		Omega_hat_4 = Omega_hat_4 + temp_Omega_hat_4
	}
	//Omega_hat_4 = min((N,T)')/T :* Omega_hat_4
		
	// Omega Hat 5th Term
	Omega_hat_5 = J(cols(X),cols(X),0)
	for( j = 1 ; j <= m-1 ; j++ ){
		temp_Omega_hat_5 = J(cols(X),cols(X),0)
		w = 1 - (j/m)
		for( utdx = (j+1) ; utdx <= T ; utdx++ ){
			tdx = uniq_t[utdx,1]
			tdxminus = uniq_t[utdx-j,1]
			for( jdx = 1 ; jdx <= cols(X) ; jdx++ ){
				for( kdx = 1 ; kdx <= cols(X) ; kdx++ ){
					temp_Omega_hat_5[jdx,kdx] = temp_Omega_hat_5[jdx,kdx] + sum((X[select((1..rows(X))',t:==tdxminus),jdx]:*Uhat[select((1..rows(X))',t:==tdxminus),1]) :* (X[select((1..rows(X))',t:==tdx),kdx]:*Uhat[select((1..rows(X))',t:==tdx),1]))
		}}} 
		temp_Omega_hat_5 = w * temp_Omega_hat_5
		Omega_hat_5 = Omega_hat_5 + temp_Omega_hat_5
	}
	
	
	// Add the 5 Terms to Get Omega Hat
	Omega_hat_4 = Omega_hat_4 + Omega_hat_4'
	Omega_hat_5 = Omega_hat_5 + Omega_hat_5'
	Omega_hat = Omega_hat_1 :+ Omega_hat_2 :- Omega_hat_3 :+ Omega_hat_4 :- Omega_hat_5
	
	
	// Q Hat
	Q_hat = X' * X
	
	//h(b) function
	bb = m/T 
	hb = 1 - bb + 1/3 * bb^2
	
	// Variance
	VCHS = luinv(Q_hat) * Omega_hat * luinv(Q_hat)
	VBCCHS = VCHS :/ hb
	VDKA = luinv(Q_hat) * (Omega_hat_1 :+ (Omega_hat_2 :+ Omega_hat_4) :/ hb) * luinv(Q_hat)
	VCi = luinv(Q_hat) * Omega_hat_1 * luinv(Q_hat)
	VCt  = luinv(Q_hat) * Omega_hat_2 * luinv(Q_hat)
	VEHW = luinv(Q_hat) * Omega_hat_3 * luinv(Q_hat)
	VDK  = luinv(Q_hat) * (Omega_hat_2 :+ Omega_hat_4) * luinv(Q_hat)
	
	//Variance type
	if (ty == 0 ){
		V= VCHS
	}
	if (ty == 1 ){
		V= VBCCHS
	}
	if (ty == 2 ){
		V= VDKA
	}
	if (ty == 3 ){
		V= VBCCHS
	}
	if (ty == 4 ){
		V= VDKA
	}
	if (ty == -1 ){
		V = (N*T)/(N*T-cols(X)) * VCi
	}
	if (ty == -2 ){
		V = (N*T)/(N*T-cols(X)) * VCt
	}
	if (ty == -3 ){
		V = (N*T)/(N*T-cols(X)) * VEHW
	}
	if (ty == -4 ){
		V = VDK
	}
	//cv output
	if (ty >= 3){
	bhat = mhat / T
	c = N / T
	hbhat = 1-bhat+1/3*bhat^2
	Lambda_ahat = 1/(N^(0.5) * T) * cholesky(Omega_hat_1)
	Lambda_ghat = 1/(N * T^(0.5)) * cholesky( (Omega_hat_2 :+ Omega_hat_4) / hbhat)
	
		t_CHS_hat = J(rep, cols(X), 0)
		dt = 1/bm
		
		R = J(cols(X), cols(X), 0)
		for (j=1; j<=cols(X); j++){
			R[j,j] = 1
		}
		
		for (k=1; k<=rep; k++){
			W = J(bm , cols(X), 0)
			
			for (j=1;j<=cols(X); j++){
				W[1,j] = invnormal(uniform(1,1))*sqrt(dt)
				for (s=2; s<=bm; s++) {
					W[s,j] = W[s-1,j] + invnormal(uniform(1,1))*sqrt(dt)
				}
			}
			
			P1 = J(cols(X),cols(X),0)
			for (r=1; r<=bm; r++) {
				P1 = P1 + 2/bb * (W[r,] - r/bm * W[bm,])' * (W[r,] - r/bm * W[bm,]) * 1/bm
			}
				
			P2 = J(cols(X),cols(X),0)
			for (r=1; r<=(1-bb)*bm; r++) {
				P2 = P2 + 1/bb * ((W[r,] - r/bm * W[bm,] )' * (W[(r/bm + bb)*bm,] - (r/bm + bb) * W[bm,] ) + (W[(r/bm + bb)*bm,] - (r/bm + bb) * W[bm,] )' * (W[r,] - r/bm * W[bm,] )) * 1/bm
			}
			P = P1 - P2
			z = J(cols(X), 1, invnormal(uniform(1,1)))
			VV = c * Lambda_ghat * P * Lambda_ghat' + hb * Lambda_ahat * Lambda_ahat'
			H = Lambda_ahat * z + c^0.5 * Lambda_ghat * W[bm,]'
			
			for (j=1;j<=cols(X); j++){	
				denominator_CHS = (R[j,] * luinv(Q_hat) * VV * luinv(Q_hat) * R[j,]')^0.5
				numerator_CHS = R[j,] * luinv(Q_hat) * H
				t_CHS_hat[k,j] = numerator_CHS / denominator_CHS	
			}
		}
		
		
		qbottom = floor(le/2*rep)
		qtop = floor((1-le/2)*rep)
	
	
		t_CHS_hat_1 = t_CHS_hat[.,whichv]
		t_CHS_hat_sort = sort(t_CHS_hat_1, 1)
		t_CHS_hat_bottom = t_CHS_hat_sort[qbottom]
		t_CHS_hat_top = t_CHS_hat_sort[qtop]
		
		t_DKA_hat_bottom = t_CHS_hat_bottom * (hb)^0.5
		t_DKA_hat_top = t_CHS_hat_top * (hb)^0.5
	
	
		cvBCCHS = (abs(t_DKA_hat_bottom) + t_DKA_hat_top)/2
		cvDKA  = cvBCCHS
	
	}
	
	if (ty < 3){
		cvBCCHS = .
		cvDKA  = .
	}
	
//////////////////////////////////////////////////////////////////////////////// 
// Output
    st_matrix(bname, beta')
    st_matrix(Vname, V)
    st_numscalar(nname, N)
	st_numscalar(tname, T)
	st_numscalar(ntname, N*T)
	st_numscalar(mname, mhat)
	st_numscalar(cvBCCHSname, cvBCCHS)
	st_numscalar(cvDKAname, cvDKA)	
}
end
////////////////////////////////////////////////////////////////////////////////

