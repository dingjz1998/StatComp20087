#include <Rcpp.h>
using namespace Rcpp;

//' @title A function to generate Markov chains using Rcpp
//' @description A function to generate Markov chains using Rcpp given the variance of the normal distribution, the initial value and the length of the chain.
//' @param sigma The variance of the normal distribution
//' @param x0 The initial value
//' @param N The length of the chain
//' @return a Markov chain with length \code{N}
//' @examples 
//' \dontrun{
//' sigma=1
//' x1=0
//' N=2000
//' X=LARMC(sigma,x0,N)
//' plot(x)
//' }
//' @export
// [[Rcpp::export]]
NumericVector LARMC(double sigma, double x0, int N) {
  NumericVector x(N);//define a function to save the chain and acp
  x[1]=x0;
  NumericVector u(N);
  u=runif(N); 
  //int acp=0 //acceptance rate
  double  y=0;
  for (int i=0;i<N;i++) { 
    y=rnorm(1, x[i-1], sigma)[0];
    if (u[i] <= (0.5*exp(-abs(y)) / (0.5*exp(-abs(x[i-1])))))
    {x[i]=y; }
    // acp=acp+1; 
    else { 
      x[i]=x[i-1];
    }
  }
  return(x);
}
