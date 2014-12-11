#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void lab(SEXP lac, SEXP lbc, SEXP SL, SEXP TA, SEXP lambda, SEXP kapp, SEXP gamm, SEXP delta, SEXP nml, int no, int n){
  using namespace Rcpp;

  //Inputs
  NumericVector SLc(SL); // Already SLc- SLminc
  NumericVector TAc(TA);
  IntegerVector nmlc(nml); // Already nmlc -1

  NumericVector lambdac(lambda);
  NumericVector kappc(kapp);
  NumericMatrix gammc(gamm);
  NumericVector deltac(delta);

  Environment CircStats("package:CircStats");
  Function dvm = CircStats["dvm"];

  NumericMatrix lalpha(lac);
  NumericMatrix lbeta(lbc);

  NumericVector foo(2);
  double sumfoo;
  double lscale;
  NumericVector foot(2);

  NumericMatrix obsProb(n,2);
  std::fill(obsProb.begin(), obsProb.end(), 1);
  NumericVector tmp(no);

  tmp = dvm(TAc,0,0);
  tmp = dexp(SLc,lambdac[0]) * tmp;
  for (int i = 0; i < no; i++){
    obsProb( nmlc[i], 0) = tmp[i];
  }
  tmp = dvm(TAc,0,kappc);
  tmp = dexp(SLc,lambdac[1]) * tmp;
  for (int i = 0; i < no; i++){
    obsProb( nmlc[i], 1) = tmp[i];
  }

  foo = deltac * obsProb( 0, _);
  sumfoo = sum(foo);
  lscale = log(sumfoo);
  foo = foo/sumfoo;
  
  lalpha( _, 0) = log(foo) + lscale;
  
  for (int i = 1; i < n; i++){
    foot[0] = foo[0] * gammc( 0, 0) + foo[1] * gammc( 1, 0);
    foot[1] = foo[0] * gammc( 0, 1) + foo[1] * gammc( 1, 1);
    foo = foot * obsProb( i, _);
    sumfoo = sum(foo);
    lscale = lscale + log(sumfoo);
    foo = foo/sumfoo;
    lalpha( _, i) = log(foo) + lscale;
  }

  foo = rep(0.5,2);
  lscale = log(2);
  for (int i = n-2; i >= 0; i--){
    foo = obsProb( i+1, _) * foo;
    foot[0] = foo[0] * gammc( 0, 0) + foo[1] * gammc( 0, 1);
    foot[1] = foo[0] * gammc( 1, 0) + foo[1] * gammc( 1, 1); 
    lbeta( _, i) = log(foot) + lscale;
    sumfoo = sum(foot);
    foo = foot/sumfoo;
    lscale = lscale + log(sumfoo);
  }

}

// [[Rcpp::export]]
SEXP EMHMM(SEXP SL, SEXP TA, SEXP ml, SEXP nml, SEXP SLmin, SEXP lambda, SEXP gamm, SEXP delta, SEXP kapp,  SEXP maxiter, SEXP tol, SEXP op){

  // Starting values for parameters to estimate
  NumericVector deltac(clone(delta));
  NumericMatrix gammc(clone(gamm));
  NumericVector lambdac(clone(lambda));
  NumericVector kappc(clone(kapp));

  // Data
  IntegerVector nmlc(clone(nml));
  IntegerVector mlc(ml);
  NumericVector SLc(clone(SL));
  NumericVector TAc(TA);
  double SLminc = as<double>(SLmin);

  // Functions from package
  Environment CircStats("package:CircStats");
  Function dvm = CircStats["dvm"];
  Function opc(op);

  // Parameters that specify the precision of the EM algorithm
  double tolc = as<double>(tol);
  int mit = as<int>(maxiter);

  // Length of the time serie, no = data, nc = true
  int no = SLc.size();
  int nc =  no + sum(mlc-1);  

  // Exponential applied to the difference in step length from minimum SL
  SLc = SLc - SLminc;
  // Indexing of C++ from 0
  nmlc = nmlc - 1;

  // Creating the matrix for the forward and backward probabilities
  NumericMatrix lac(2,nc);
  NumericMatrix lbc(2,nc);
  // Creating matrix for the log probabilities of observation
  NumericMatrix lprobObsc(nc,2);

  // Creating the variable for the log-likelihood
  double llkc;

  // Set of temporary variables used to get at estimated parameters
  NumericVector term2(1);

  NumericMatrix gammn(2,2);
  NumericVector lambdan(2);
  NumericVector deltan(2);
  NumericVector kappan(1);

  NumericVector gg(2); // Because matrices are difficult to handle

  NumericVector tmp1(nc-1);
  NumericVector tmp2(nc-1);
  NumericVector tmp3(no);
  NumericVector tmp4(no);

  List res = List::create(Named("lambda"), Named("gamma"), Named("delta"), Named("kapp"), Named("mllk"));

  double crit;

  for (int ii =1; ii <= mit; ii++){

    lab(lac, lbc, SLc, TAc, lambdac, kappc, gammc, deltac, nmlc, no, nc);

    llkc = fmax(lac(0,nc-1),lac(1,nc-1));
    
    if(R_IsNaN(lac( 0, nc-1)) || R_IsNaN(lac( 1, nc-1))){
        Rf_warning("Underflow problem cannot calculate the log likelihood.");
        res["lambda"] = NA_REAL;
        res["gamma"] = NA_REAL;
        res["kapp"] = NA_REAL;
        res["delta"] = NA_REAL;
        res["mllk"] = NA_REAL;
        break;
    }

    llkc = llkc + log(exp(lac(0,nc-1) - llkc) + exp(lac(1,nc-1) - llkc));


    tmp3 = dvm(TAc,0,0);
    tmp3 = dexp(SLc,lambdac[0],1) + log(tmp3);
    tmp4 = dvm(TAc,0,kappc[0]);
    tmp4 = dexp(SLc,lambdac[1],1) + log(tmp4);  
    for (int j = 0; j < no; j++){
      lprobObsc( nmlc[j], 0) = tmp3[j];
      lprobObsc( nmlc[j], 1) = tmp4[j];
    } 

    for (int j = 0; j < 2; j++){
      for (int k =0; k < 2; k++){
        tmp1 = lac(Range(j,j),Range(0,(nc-2)));
        tmp2 = lprobObsc(Range(1,(nc-1)),Range(k,k));
        tmp1 = tmp1 + tmp2;
        tmp2 = lbc(Range(k,k),Range(1,(nc-1)));
        tmp1 = tmp1 + tmp2 - llkc;
        tmp1 = exp(tmp1);
        gammn(j,k) = gammc(j,k) * sum(tmp1);
      }
      for (int k = 0; k < no; k++){
        tmp3[k] = lac( j, nmlc[k]) + lbc( j, nmlc[k]);
      }
      tmp3 = tmp3 - llkc;
      tmp3 = exp(tmp3);
      tmp4 = tmp3 * SLc;
      lambdan[j] = sum(tmp3) / sum(tmp4);
    }
    gammn( 0, _) = gammn( 0, _) / sum(gammn( 0, _));
    gammn( 1, _) = gammn( 1, _) / sum(gammn( 1, _));

    deltan = lac( _, 0) + lbc( _, 0) - llkc;
    deltan = exp(deltan);
    deltan = deltan / sum(deltan);

    // Estimating the CDLL of kappa_T
    for (int j = 0; j < no; j++){
      tmp3[j] = lac( 1, nmlc[j]) + lbc( 1, nmlc[j]) - llkc;
    }
    tmp3 = exp(tmp3);
    tmp4 = tmp3 * cos(TAc);

    term2 = sum(tmp4) / sum(tmp3);

    kappan = opc(term2);
  
    crit = sum(abs(lambdac - lambdan)) + sum(abs(deltac - deltan)) + sum(abs(kappc - kappan));
    gg = gammc( _, 0) - gammn( _, 0);
    gg = abs(gg);
    crit = crit + sum(gg);
    gg = gammc( _, 1) - gammn( _, 1);
    gg = abs(gg);
    crit = crit + sum(gg);
    if (crit < tolc){
      res["lambda"] = lambdac;
      res["gamma"] = gammc;
      res["kapp"] = kappc;
      res["delta"] = deltac;
      res["mllk"] = -llkc; 
      break;
    }
  
    if(is_true(any(lambdan == R_PosInf))){
      Rf_warning("Potential unbounded likelihood problem. One of the lambda has an infinite value.");
      res["lambda"] = NA_REAL;
      res["gamma"] = NA_REAL;
      res["kapp"] = NA_REAL;
      res["delta"] = NA_REAL;
      res["mllk"] = NA_REAL;
      break;
    }
    
    lambdac = clone(lambdan);
    gammc = clone(gammn);
    deltac = clone(deltan);
    kappc = clone(kappan);

    if (ii == mit){
      Rf_warning("No convergence after the specified maximum number of iterations");
      res["lambda"] = NA_REAL;
      res["gamma"] = NA_REAL;
      res["kapp"] = NA_REAL;
      res["delta"] = NA_REAL;
      res["mllk"] = NA_REAL;
      break;
    }
  }
  return(res);

}
