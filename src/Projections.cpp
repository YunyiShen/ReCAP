// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <climits>
#include <complex>
using namespace Rcpp;

///get Leslie matrix from survival etc.
// [[Rcpp::export]]
arma::mat getLeslie(const arma::mat& Surv, const arma::mat& Fec, const double& SRB){
	int nage_female = Fec.n_elem;
	int nage_male = Surv.n_elem - nage_female;
	arma::mat Tr(nage_female + nage_male , nage_male + nage_female);
	Tr.zeros();
	Tr.row(0).cols(0,nage_female-1) = SRB * Fec.t();
	Tr.row(nage_female).cols(0,nage_female-1) = (1 - SRB) * Fec.t();
	Tr.submat(1,0,nage_female-1,nage_female-2).diag() = Surv.rows(0,nage_female-2).t();
	Tr(nage_female-1,nage_female-1) = Surv(nage_female-1,0);
	if(nage_male>1) Tr.submat(nage_female + 1,nage_female,nage_male + nage_female-1,nage_male + nage_female-2).diag() = Surv.rows(nage_female,nage_male + nage_female-2).t();
	Tr(nage_male+nage_female-1,nage_male+nage_female-1) = Surv(nage_male+nage_female-1,0);
	return(Tr);
}

///Calculate the density dependency
arma::mat DD(const bool& global, const arma::mat& Xn,const arma::mat & aK0, const arma::mat& midP, const bool& null){
  //E0 = E0/sum(E0);// This was done in main projector
  arma::mat D;
  if(global){
    arma::mat den = ( 1 + (aK0) * (sum(Xn)) )/(1+aK0 * midP);
    D = (1-null)*den + null;
  }
  else{
	arma::mat den = ( 1 + (aK0) % ((Xn)) )/(1+aK0 % midP);
    D = (1-null)*den + null;
  }
  return(D);
}

///Helper function for a single year projection, inner function, export for test.
arma::mat ProjectHarvest_helperCpp(const arma::mat& data_n,const arma::mat& Surv, const arma::mat& Fec,const double& SRB,const arma::mat& H_n, const arma::mat& H_np1,bool global, const List& aK0,const bool & null){
	arma::mat X_n1 = (1-H_n) % (data_n/H_n);
	arma::mat D_bir = DD(global, X_n1, aK0[0], aK0[2] ,null);
	arma::mat D_dea = DD(global, X_n1, aK0[1], aK0[2] ,null);
	return(H_np1 % (getLeslie(Surv % D_dea, Fec % D_bir, SRB)*X_n1));

}

///main projection function
//[[Rcpp::export]]
arma::mat ProjectHarvestCpp(const arma::mat& Surv,const arma::mat& Harvpar,const arma::mat& Fec, const arma::mat& SRB, const List& aK0, const bool& global, const bool& null, const arma::mat& bl ,const int& period, const IntegerVector& nage){
	arma::mat Harvest(sum(nage),period+1);
	Harvest.col(0) = bl;
	//E0 = E0/(sum(E0));// need to check whether there is one in R call rather than here.
	for(int i = 1; i<period + 1; i++){
		Harvest.col(i) = ProjectHarvest_helperCpp(Harvest.col(i-1),Surv.col(i-1),Fec.col(i-1),(SRB(0,i-1)), Harvpar.col(i-1),Harvpar.col(i),global, aK0,null);
	}
	return(Harvest);
}

///get Aerial count
//[[Rcpp::export]]
arma::mat getAerialCount(const arma::mat& Harv, const arma::mat& H, const arma::mat& A){
  return((sum((1/H-1) % Harv))%A);
}


///get observed lambda (growth rate from harvest)

arma::mat get_obs_LambdasH(const arma::mat& Harv, const arma::mat& H){
  arma::mat Living_total = sum((1/H-1) % Harv);//living individual at all year
  return((Living_total.cols(0,Harv.n_cols-2))/(Living_total.cols(1,Harv.n_cols-1))); // Lambdas
}

///get observed lambda (growth rate from Aerial count)
//[[Rcpp::export]]
arma::mat get_obs_LambdasA(const arma::mat& Ae, const arma::mat& Ae_det){
  arma::mat Living_total = Ae/Ae_det;//living individual at all year
  return((Living_total.cols(1,Ae.n_cols-1))/(Living_total.cols(0,Ae.n_cols-2))); // Lambdas
}

///get lambda w/, w/o harvest and maximum lambda, single year
//[[Rcpp::export]]
arma::mat get_hypo_Lambdas_helper( const arma::mat& Harv_n // harvest count
                              , const arma::mat& H_n // harvest rate at year n
                              , const arma::mat& Surv_np1 // survival of year n+1
                              , const arma::mat& Fec_np1 // Fecundity of year n+1
                              , const double& SRB_np1 // SRB of year n+1
                              , const arma::mat& H_np1){
  arma::mat res(5,1);
  arma::mat living = (1/H_n-1) % Harv_n;
  arma::mat Leslie_np1 =  getLeslie(Surv_np1,Fec_np1,SRB_np1);
  // eigen problem, possible lambda

  //max intrinsic
  res.row(0).col(0) = max(sum(Leslie_np1));

  arma::mat avg_age_str(H_n.n_rows,1);
  avg_age_str.rows(0,Fec_np1.n_rows-1) = 0.5 * 1/(Fec_np1.n_rows);
  avg_age_str.rows(Fec_np1.n_rows,Surv_np1.n_rows-1) = 0.5 * 1/(Surv_np1.n_rows-Fec_np1.n_rows);
  //even age structure
  res.row(1).col(0) = sum(Leslie_np1 % avg_age_str);// uniform age structure, assuming 1 to 1 sex ratio
  //stable intrinsic
  res.row(2).col(0) = real(max((eig_gen(Leslie_np1))));// get the intrinsic one

  // w/o harvest:

  res.row(3).col(0) = sum(Leslie_np1 * living)/sum(living);

  //min intrinsic
  res.row(4).col(0) = min(sum(Leslie_np1));// basically when population are all female fawns.

  return(res);

}

///get lambda w/, w/o harvest and maximum lambda, single year
//[[Rcpp::export]]
arma::mat get_hypo_Lambdas(const arma::mat& Harvest
                             , const arma::mat& Harvpar
                             , const arma::mat& Surv
                             , const arma::mat& Fec
                             , const arma::mat& SRB){
  int periods = Harvest.n_cols-1;
  arma::mat Lambdas(5,periods);
  for(int i=1;i<periods+1;i++){
    Lambdas.col(i-1)=get_hypo_Lambdas_helper( Harvest.col(i-1) , Harvpar.col(i-1) , Surv.col(i-1) , Fec.col(i-1) , SRB(0,i-1) , Harvpar.col(i));
  }
  return(Lambdas);
}

///Misc
//[[Rcpp::export]]
arma::mat eyes(const int& n){
  arma::mat I;
  I.eye(n,n);
  return(I);
}
