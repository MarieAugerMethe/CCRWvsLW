#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP sslta(SEXP n, SEXP a, SEXP lambda, SEXP kapp, SEXP ai2, SEXP delta){
 
  int nc = as<int>(n);
  double ac = as<double>(a);
  NumericVector kappc(kapp);
  double deltac = as<double>(delta);
  NumericVector lambdac(lambda);
  NumericVector ai2c(ai2);

  NumericVector bpc(nc+1);
  NumericMatrix coord(nc,2);
  NumericVector dist(1);
  NumericVector beh(1);

  // I rewrote the rvm function these are variable needed for it
  double U1;
  double U2;
  double U3;
  double f;
  double c;
  int i;  

  double b = 1 + pow((1 + 4 * pow(kappc[0],2)),0.5);
  b = (b - pow(2*b,0.5))/(2*kappc[0]);
  double r = (1 + pow(b,2))/(2 * b);

  RNGScope scope;


  // Initialize beh Markov chain
  bpc[0] = as<double>(rbinom(1, 1, deltac));

  for(int j = 0; j < nc; j++){
    dist = rexp(1,lambdac[bpc[j]]);
    coord(j,0) = dist[0] + ac;
    if(bpc[j] == 0){
      coord(j,1) = as<double>(runif(1,0,2*PI));
    }
    else{
      i = 0;
      while(i < 1){
        U1 = as<double>(runif(1, 0, 1));
        f = cos(PI * U1);
        f = (1 + r * f)/(r + f);
        c = kappc[0] * (r - f);
        U2 = as<double>(runif(1, 0, 1));
        if(c * (2 - c) - U2 > 0 || log(c/U2) + 1 - c >= 0){
          U3 = as<double>(runif(1, 0, 1));
          if(U3 < 0.5){
            coord(j,1) = 2 * PI - acos(f);
          }else{
            coord(j,1) = acos(f);
          }
            coord(j,1) = fmod(coord(j,1), 2*PI);
            i = 1;
          }
        }
      }

      beh = rbinom(1,1,ai2c[bpc[j]]);
      bpc[j+1] = beh[0];
    }

  return coord;
}
