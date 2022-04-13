//simulate.cpp
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double rcpp_sum(NumericVector v){
    double sum = 0;
    for(int i=0; i<v.length(); ++i){
        sum += v[i];
    }
    return(sum);
}


// [[Rcpp::export]]
NumericVector rcpp_rbinom(int n, int size, double prob) { 
    NumericVector out = rbinom(n, size, prob);
    return(out);
}


// [[Rcpp::export]]
NumericVector rcpp_rpois(int n, double lambda) { 
    NumericVector out = rpois(n, lambda);
    return(out);
}


// [[Rcpp::export]]
IntegerVector rcpp_mutate(IntegerVector gt, double mu, NumericVector biases) {
    int n_markers=gt.length();
    NumericVector indels = rbinom( n_markers, 1, mu );

    for(int i=0; i<n_markers; i++){

        // for the currrent marker, get indel direction considering the bias
        // NB: 0 means deletion, 1 means insertion
        double insertion = R::rbinom( 1, biases[i] );

        if(insertion==0){
            indels[i] = -1*indels[i];
        }
     
        // if the value at pos i is not NA, then update it
        if(IntegerVector::is_na(gt[i]) == false) {
            gt[i] = gt[i] + indels[i];           
        }
     
    }
    return(gt);
}


// [[Rcpp::export]]
IntegerMatrix rcpp_mutate_length_matrix(IntegerMatrix gt, double mu, int gens, NumericVector biases) { 
    int max_ploidy = gt.nrow(); 

    for(int j=0; j<max_ploidy; j++) {
        for(int i=0; i<gens; i++){
            gt(j,_) = rcpp_mutate(gt(j,_), mu, biases);
        }
    }
    return(gt);
}




