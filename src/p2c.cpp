#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP p2c(SEXP SL, SEXP TA){
  // Input values
  Rcpp::NumericVector distc(SL);
  Rcpp::NumericVector rtac(TA);
  int nc = distc.size();
  nc += 1;
  // Create dx and dy vector
  NumericVector dx(nc);
  NumericVector dy(nc);
  // Create coord matrix and set the initial location to 0
  Rcpp::NumericMatrix coord(nc,2);

  // Transforming relative turning angle into absolute turning angle
  NumericVector ata(nc);
  ata[0] = fmod(rtac[0],2*PI);
  for(int j = 1; j < nc; j++){  
    ata[j] = fmod(ata[j-1] + rtac[j],2*PI);
  }

  /////////////////////////////////
  //Change into cartesian coordinates

  // Calculate the displacement in each direction

  dx = distc * cos(ata);
  dy = distc * sin(ata);

  // Add the displacments
  for(int j = 1; j < nc; j++){
    coord(j,0) = coord(j-1, 0) + dx[j-1];
    coord(j,1) = coord(j-1, 1) + dy[j-1];
  }

  return coord;
}