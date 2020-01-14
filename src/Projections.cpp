// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <climits>
#include <complex>
using namespace Rcpp;

arma::mat ReLU( arma::mat matr){
    int n_elem = matr.n_elem;
    for(int i = 0; i < n_elem ; ++i){
        matr(i) = std::max(matr(i),0.0);
    }
    return(matr);
}

///get Leslie matrix from survival etc.
// [[Rcpp::export]]
arma::mat getLeslie(const arma::mat& Surv,
					const arma::mat& Fec,
					const double& SRB){
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
arma::mat DD(const bool& global,
			 const arma::mat& Xn,
			 const arma::mat & aK0,
			 const arma::mat& midP,
			 const bool& null){
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


Rcpp::List ProjectLiving_helper2Cpp(const arma::mat& Living_n1,
									const arma::mat& Surv,
									const arma::mat& Fec,
									const double& SRB,
									const bool& global,
									const List& aK0,
									const bool & null){
	arma::mat D_bir = DD(global, Living_n1, aK0[0], aK0[2] ,null);
	arma::mat D_dea = DD(global, Living_n1, aK0[1], aK0[2] ,null);
	arma::mat Fec_obs = Fec % D_bir;
	arma::mat Surv_obs = Surv % D_dea;
	arma::mat Living = (getLeslie(Surv_obs, Fec_obs, SRB)*Living_n1);
	return(Rcpp::List::create(
		Rcpp::Named("Living")=Living,//pre harvest living individuals
		Rcpp::Named("Fec") = Fec_obs,// observed fecundity w/ DD
		Rcpp::Named("Surv") = Surv_obs// observed survival w/ DD
		)
	);
}


///main projection function, w/ missing harvest
//[[Rcpp::export]]
Rcpp::List ProjectAllCpp(const arma::mat& Surv,
						 const arma::mat& Harvpar,
						 const arma::mat& Fec,
						 const arma::mat& SRB,
						 const List& aK0,
						 const bool& global,
						 const bool& null,
						 const arma::mat& bl ,
						 const int& period,
						 const IntegerVector& nage){
	arma::mat Harvest(sum(nage),period+1); // deal with the case that there is no harvest for a certain year
	arma::mat Living(sum(nage),period+1);
	arma::mat Fecobs(nage(0),period);
	arma::mat Survobs(sum(nage),period);
	arma::mat living_temp(sum(nage),1);
	List Proj_temp(3);
	Harvest.col(0) = bl;
	Living.col(0) = (1/Harvpar.col(0)) % bl; // pre harvest living individuals
	for(int i = 1; i<period + 1; i++){
			Proj_temp = ProjectLiving_helper2Cpp(Living.col(i-1)-Harvest.col(i-1), Surv.col(i-1),Fec.col(i-1),(SRB(0,i-1)),global, aK0,null);// pre-cull population of next year
			arma::mat Living_temp = Proj_temp["Living"];
			Living.col(i) = Living_temp;
			arma::mat Fec_temp = Proj_temp["Fec"];
			Fecobs.col(i-1) = Fec_temp;
			arma::mat Surv_temp = Proj_temp["Surv"];
			Survobs.col(i-1) = Surv_temp;
			Harvest.col(i) = Living.col(i) % Harvpar.col(i); // culling
	}
	//List res = List::create(Named["Harvest"] = Harvest , _["Living"] = Living);
	return(Rcpp::List::create(
	        Rcpp::Named("Harvest") = Harvest ,
	        Rcpp::Named("Living") = Living-Harvest,
		  	Rcpp::Named("Fec_obs") = Fec,
		    Rcpp::Named("Surv_obs") = Surv
	));//
}

///get Aerial count
//[[Rcpp::export]]
arma::mat getAerialCountPost(const List& Proj,
							 const arma::mat& obsMat,
							 const arma::mat& A){
  arma::mat living = Proj["Living"];
  return(((obsMat.t() * (living%A))));
}
///get Aerial count
//[[Rcpp::export]]
arma::mat getAerialCountPre(const List& Proj,
							const arma::mat& obsMat,
							const arma::mat& A){
    arma::mat living = Proj["Living"];
    arma::mat harv = Proj["Harvest"];
    return((obsMat.t() * ((living+harv)%A)));//
}

//[[Rcpp::export]]
arma::mat getobsVitals(const arma::mat& vital,
					   const arma::mat& living,
					   const arma::mat& obsMat){
	arma::mat previous_living = living.cols(0,living.n_cols-2);
    arma::mat res(obsMat.n_cols,living.n_cols-1);
    for(int i = 0 ; i<living.n_cols-1;++i){
        res.col(i) = (obsMat.t() * (vital.col(i) % previous_living.col(i)))/
            (obsMat.t() * previous_living.col(i));
    }
	return(res); // matrix form for weighted average of vital rate,
// classMat has row as age classes, col as classes, e.g. 5 age and 3 classes, first two are two classes and last 3 are one, classMat will be
	//[1,0,0;0,1,0;0,0,1;0,0,1;0,0,1]
}

///get observed lambda (growth rate from harvest)

arma::mat get_obs_LambdasH(const arma::mat& Harv,
						   const arma::mat& H){
  arma::mat Living_total = sum((1/H-1) % Harv);//living individual at all year
  return((Living_total.cols(0,Harv.n_cols-2))/(Living_total.cols(1,Harv.n_cols-1))); // Lambdas
}

///get observed lambda (growth rate from Aerial count)
//[[Rcpp::export]]
arma::mat get_obs_LambdasA(const arma::mat& Living_total){
  return((Living_total.cols(1,Living_total.n_cols-1))/(Living_total.cols(0,Living_total.n_cols-2))); // Lambdas
}

///get lambda w/, w/o harvest and maximum lambda, single year
//[[Rcpp::export]]
arma::mat get_hypo_Lambdas_helper( const arma::mat& Harv_n // harvest count
							  , const arma::mat& living // post cull living individual
                              , const arma::mat& H_n // harvest rate at year n
                              , const arma::mat& Surv_np1 // survival of year n+1
                              , const arma::mat& Fec_np1 // Fecundity of year n+1
                              , const double& SRB_np1 // SRB of year n+1
                              , const arma::mat& H_np1){
  arma::mat res(5,1);
  arma::mat Leslie_np1 =  getLeslie(Surv_np1,Fec_np1,SRB_np1);
  // eigen problem, possible lambda

  //max intrinsic
  res.row(0).col(0) = max(sum(Leslie_np1));

  arma::mat avg_age_str(H_n.n_rows,1);
  avg_age_str.ones();
  avg_age_str.rows(0,Fec_np1.n_rows-1) = (0.5 * 1/(Fec_np1.n_rows)) * avg_age_str.rows(0,Fec_np1.n_rows-1);
  avg_age_str.rows(Fec_np1.n_rows,Surv_np1.n_rows-1) = (0.5 * 1/(Surv_np1.n_rows-Fec_np1.n_rows)) * avg_age_str.rows(Fec_np1.n_rows,Surv_np1.n_rows-1);
  //even age structure
  res.row(1).col(0) = sum(Leslie_np1 * avg_age_str);// uniform age structure, assuming 1 to 1 sex ratio
  //stable intrinsic
  res.row(2).col(0) = real(max((eig_gen(Leslie_np1))));// get the intrinsic one

  // w/o harvest:

  res.row(3).col(0) = sum(Leslie_np1 * living)/sum(living);

  //min intrinsic
  res.row(4).col(0) = min(sum(Leslie_np1));// basically when population are all female fawns.

  return(res);
}



// Harvest with certain weight and quota
arma::mat get_hypoharv_quota_helper( const arma::mat& living_n // post cull living individual of year n
                              , const arma::mat& Surv_np1 // survival of year n+1
                              , const arma::mat& Fec_np1 // Fecundity of year n+1
                              , const double& SRB_np1 // SRB of years
							  , const bool& global
							  ,	const List& aK0
							  , const bool & null
							  , const arma::mat& Harv_np1 // # harvest of year n+1
							  , const arma::mat & hypoharv_np1 // harvest weight, to test different scheme, e.g. solely adult female can be 0,0,1,1,1,1,1,1,0,0,0
                              ){
	arma::mat D_bir = DD(global, living_n, aK0[0], aK0[2] ,null);
	arma::mat D_dea = DD(global, living_n, aK0[1], aK0[2] ,null);
	arma::mat Fec_obs = Fec_np1 % D_bir;
	arma::mat Surv_obs = Surv_np1 % D_dea;
	arma::mat Living = (getLeslie(Surv_obs, Fec_obs, SRB_np1)*living_n); //just let the population increase

	// now we need to harvest:
	arma::mat harv_np1 =  ((Living % hypoharv_np1) * (1/(sum(Living % hypoharv_np1)))) * sum(Harv_np1) ;

	return( ReLU( Living-harv_np1));

}

// Harvest quota scheme
//[[Rcpp::export]]

arma::mat get_hypo_harvest_quotaCpp(const arma::mat& bl // baseline before harvest
						  , const arma::mat& Harv
						  , const arma::mat& Surv
						  , const arma::mat& Fec
						  , const arma::mat& SRB
						  , const arma::mat& harv_weight
						  , const bool& global
						  ,	const List& aK0
						  , const bool & null
						  ){
	arma::mat Living = 0 * Harv;
    Living.col(0) = bl - ((bl % harv_weight.col(0)) * (1/(sum(bl % harv_weight.col(0))))) * sum(Harv.col(0));
	Living.col(0) = ReLU(Living.col(0)) ;//first year, baseline

	int period = Harv.n_cols - 1;

	for(int i=1; i<period + 1; ++i){
		Living.col(i) = get_hypoharv_quota_helper(Living.col(i-1) , Surv.col(i-1) , Fec.col(i-1) , SRB(0,i-1) , global , aK0 , null , Harv.col(i) , harv_weight.col(i));
	}

	return(Living);
}

// Harvest with certain weight and portion
//[[Rcpp::export]]
arma::mat get_hypo_harvest_portionCpp(const arma::mat& bl // baseline before harvest

                                        , const arma::mat& Harv_rate // this should be harvest rate in whole population level
                                        , const arma::mat& Surv
                                        , const arma::mat& Fec
                                        , const arma::mat& SRB
                                        , const arma::mat& harv_weight
                                        , const bool& global
                                        ,	const List& aK0
                                        , const bool & null
                                        , const int &period
                                        , const IntegerVector& nage
){
    arma::mat Harv = sum( bl) * Harv_rate(0,0);
    arma::mat Living(sum(nage),period+1);
    Living.col(0) = bl - ((bl % harv_weight.col(0)) * (1/(sum(bl % harv_weight.col(0))))) * Harv;
    Living.col(0) = ReLU(Living.col(0)) ;//first year, baseline

    for(int i=1; i<period + 1; ++i){
        Harv = sum(Living.col(i-1)) * Harv_rate(i,0);
        Living.col(i) = get_hypoharv_quota_helper(Living.col(i-1) , Surv.col(i-1) , Fec.col(i-1) , SRB(0,i-1) , global , aK0 , null , Harv , harv_weight.col(i));
    }

    return(Living);
}






///get lambda w/, w/o harvest and maximum lambda, single year
//[[Rcpp::export]]
arma::mat get_hypo_Lambdas(const arma::mat& Harvest
							 , const arma::mat& Living
                             , const arma::mat& Harvpar
                             , const arma::mat& Surv
                             , const arma::mat& Fec
                             , const arma::mat& SRB){
  int periods = Harvest.n_cols-1;
  arma::mat Lambdas(5,periods);
  for(int i=1;i<periods+1;i++){
    Lambdas.col(i-1)=get_hypo_Lambdas_helper( Harvest.col(i-1) , Living.col(i-1),Harvpar.col(i-1) , Surv.col(i-1) , Fec.col(i-1) , SRB(0,i-1) , Harvpar.col(i));
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



