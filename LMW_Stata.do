
/* Stochastic Dominance Test of Linton-Maasoumi-Whang(2005) */

cap program drop lmwtest
program lmwtest ,eclass 
	syntax varlist(min=2 max=2) [,b(integer 10) grid(integer 10) s(integer 10)]
	ereturn clear
	mata : lmwtest("`varlist'", B = `b', grid_num = `grid', s=`s')  
	ereturn scalar test_stat = lmw
	ereturn scalar pvalue = pv
end

mata
mata clear
real matrix lmwtest(varlist, B, grid_num, s){	
	st_view(data=.,.,varlist)
	x1 = data[,1]
	x2 = data[,2]
	n1 = sum(rowmissing(data[,1]):==0)
	n2 = sum(rowmissing(data[,2]):==0)
	x_grid = rangen(min(x1\x2),max(x1\x2),grid_num)

	lmw = stat(x1,x2,x_grid,s)
	pval = mean(subsampling(x1,x2,x_grid,1,10,10) :> lmw)

	st_numscalar("lmw", lmw)
	st_numscalar("pv", pval)
	//st_matrix("pv", pval)
		
	/* Print Result */
	printf("{txt} LMW Test for Stochastic Dominance \n")
	printf("{txt} * H0 : x1 FSD x2 \n")
	printf("{txt} * Number of grid \t= %5.0g \n",grid_num)
	printf("{txt} * Number of bootstrap \t= %5.0g\n",B)
	printf("{txt} * Dominance order \t= %5.0g\n\n",s)
	
	printf("{txt}{space 25}{c |} {space 3} LMW {space 5} p-value {space 5}\t\n")
	printf("{hline 25}{c +}{hline 26}\t\n")
	printf("{txt} %20s \t {c |} {space 2} %5.4f {space 3} %5.4f \t\n","Multiplier method1",lmw,pval)
	//return(lmw,pval)
	
}
real matrix oper(X,z,s){
	Z = J(rows(X),1,z')
	op = J(rows(X),rows(z),0)
		
	for (j=1; j<=rows(z); j++){
		op[,j] = (X :<= Z[,j]):*(Z[,j]-X):^(s-1):/ factorial(s-1)
	}
	return(op)
}
real matrix ecdf(X,z,s){
	return(mean(oper(X,z,s)))
}
real matrix stat(x1,x2,z,s){
	n1 = rows(x1)
	n2 = rows(x2)
	stat = sqrt(n1*n2/(n1+n2))*max(ecdf(x1,z,s)-ecdf(x2,z,s))
	return(stat)
}


real matrix subsample(x1,x2,z,s,b){
	n1 = rows(x1)
	n2 = rows(x2)
	n = rows(x1)
	nsub = n-b+1
	subindex = J(b,nsub,0)
	for (j=1; j<=nsub; j++){
		subindex[j,] = rangen(0,nsub-1,nsub)' :+ j
	}
	stat = J(k,1,0)
	for (k=1; k<=nsub; k++){
		sqrt(b)*max(ecdf(x1[subindex[,k]],z,s):-ecdf(x2[subindex[,k]],z,s))
	}
	return(k)
}
real matrix subsampling(x1,x2,z,s,b1,b2){
	n1 = rows(x1)
	n2 = rows(x2)
	n = n1*n2/(n1+n2)
	lambda = n2/(n1+n2)
	b = lambda * b1
	nsub = min(n1-b1+1\n2-b2+1)
	
	subindex1 = J(b1,nsub,0)
	for (j=1; j<=b1; j++){
		subindex1[j,] = rangen(0,nsub-1,nsub)' :+ j 
	}
	subindex2 = J(b2,nsub,0)
	for (j=1; j<=b2; j++){
		subindex2[j,] = rangen(0,nsub-1,nsub)' :+ j 
	}
	stat = J(nsub,1,0)
	for (k=1; k<=nsub; k++){
		stat[k,] = sqrt(b)*max(ecdf(x1[subindex1[,k]],z,1):-ecdf(x2[subindex2[,k]],z,1))
	}
	return(stat)
}	
end


/******************************************************************************/
* Test
* H0 : x1 FSD x2 // H1 : x1 not FSD x2
* x1 ~ N(-0.5,1), x2 ~ N(0,1) : x2 FSD,SSD x1
* True : reject H0 for s = 1,2

*import excel "C:\data_temp\ex.xlsx",sheet("Sheet1") first clear
*Test
clear
global n 500
set obs $n

gen x1 = rnormal(0,1)
gen x2 = rnormal(0,1)

lmwtest x1 x2,b(10) grid(50) s(1)




/* Simulation
capture program drop lmwsim
program lmwsim,rclass
clear
global n 500
set obs $n

gen x1 = rnormal(0,1)
gen x2 = rnormal(0.5,1)

lmwtest x1 x2,b(10) grid(50) s(1)
return scalar pv = e(pvalue) < 0.05
end

lmwsim
simulate a = r(pv), reps(100) seed(123447) nodots : lmwsim
su
