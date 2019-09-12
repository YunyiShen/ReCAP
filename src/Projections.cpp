// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <climits>
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
	Tr.submat(nage_female + 1,nage_female,nage_male + nage_female-1,nage_male + nage_female-2).diag() = Surv.rows(nage_female,nage_male + nage_female-2).t();
	Tr(nage_male+nage_female-1,nage_male+nage_female-1) = Surv(nage_male+nage_female-1,0);
	return(Tr);
}

///Calculate the density dependency
arma::mat DD(const bool& global, const arma::mat& Xn,const arma::mat & aK0, const arma::mat& intc, const bool& null){
  //E0 = E0/sum(E0);// This was done in main projector
  arma::mat D;
  if(global){
    arma::mat D = intc + (aK0 * (1-null))*sum(Xn);

  }
  else{
	arma::mat D = intc + (aK0 * (1-null)) % ((Xn)) ;

  }
  return(D);
}

///Helper function for a single year projection, inner function, export for test.
arma::mat ProjectHarvest_helperCpp(const arma::mat& data_n,const arma::mat& Surv, const arma::mat& Fec,const double& SRB,const arma::mat& H_n, const arma::mat& H_np1){
	arma::mat X_n1 = (1-H_n) % (data_n/H_n);
	return(H_np1 % (getLeslie(Surv, Fec, SRB)*X_n1));

}

arma::mat ProjectDDHarvest_helperCpp(const arma::mat& data_n,const double& SRB,const arma::mat& H_n, const arma::mat& H_np1,const List& aK0, const bool & global, const bool & null){
	arma::mat X_n1 = (1-H_n) % (data_n/H_n);
	arma::mat Fec = exp(DD(global,X_n1,aK0[0],aK0[1],null));
	arma::mat logit_Surv = DD(global,X_n1,aK0[2],aK0[3],null);
	arma::mat Surv = exp(logit_Surv)/(1+exp(logit_Surv));
	return(H_np1 % (getLeslie(Surv, Fec, SRB)*X_n1));

}


///main projection function
//[[Rcpp::export]]
arma::mat ProjectHarvestCpp(const arma::mat& Surv,const arma::mat& Harvpar,const arma::mat& Fec, const arma::mat& SRB,const bool& lm_vital, const List& aK0, const bool& global, const bool& null, const arma::mat& bl ,const int& period, const IntegerVector& nage){
	arma::mat Harvest(sum(nage),period+1);
	Harvest.col(0) = bl;
	if(!lm_vital){
		for(int i = 1; i<period + 1; i++){
			Harvest.col(i) = ProjectHarvest_helperCpp(Harvest.col(i-1),Surv.col(i-1),Fec.col(i-1),(SRB(0,i-1)), Harvpar.col(i-1),Harvpar.col(i));
		}
	}
	else{
		for(int i = 1; i<period + 1; i++){
			Harvest.col(i) = ProjectDDHarvest_helperCpp(Harvest.col(i-1),SRB(0,i-1), Harvpar.col(i-1),Harvpar.col(i),aK0,global,null);//aK0 here should assumed to be full version, done in R
		}
	}
	return(Harvest);
}

///get Aerial count
//[[Rcpp::export]]
arma::mat getAerialCount(const arma::mat& Harv, const arma::mat& H, const arma::mat& A){
  return((sum((1/H-1) % Harv))%A);
}

/// get linear survival and fecundity
//[[Rcpp::export]]
List getDDVitalRate(const arma::mat& Harv, const arma::mat& H, const List & aK0,const bool & global, const bool & null, const IntegerVector& nage){
	arma::mat X = (1/H-1) % Harv; // living individual before harvest
	X = X.cols(0,X.n_cols -2); // we do not need last year
    List Vital(2);
    arma::mat log_Fec(nage(0),X.n_cols);
    arma::mat logit_Surv(sum(nage),X.n_cols);
	for(int i=0; i<X.n_cols;i++){
    log_Fec.col(i) = DD(global,X.col(i),aK0[0],aK0[1],null);
    logit_Surv.col(i) = DD(global,X.col(i),aK0[2],aK0[3],null);
	}
	Vital[1] = exp(log_Fec);
	Vital[2] = exp(logit_Surv)/(1+exp(logit_Surv));

	return(Vital);
}


///Misc
//[[Rcpp::export]]
arma::mat eyes(const int& n){
  arma::mat I;
  I.eye(n,n);
  return(I);
}


