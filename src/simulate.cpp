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
NumericVector rcpp_mutate(NumericVector gt, double mu) { 
    int n_markers=gt.length();
    
    // decide which markers have an indel 
    NumericVector indels = rcpp_rbinom( n_markers, 1, mu );

    // assuming each marker had an indel (assume marker equally likely to be insertion or deletion)
    NumericVector change_in_length = rcpp_rbinom( n_markers, 1, 0.5); 

    // change in length is either -1 (deletion) or +1 (insertion)
    change_in_length = (change_in_length*2)-1;

    // markers without indels have 0 change in length 
    change_in_length = change_in_length*indels;

    for(int i=0; i<n_markers; i++){
        // if the value at pos i is not NA, then update it
        if(IntegerVector::is_na(gt[i]) == false) {
            gt[i] = gt[i] + change_in_length[i];           
        }
    }
    return(gt);
}


// [[Rcpp::export]]
NumericMatrix rcpp_mutate_length_matrix(NumericMatrix gt, double mu, int gens) { 
    int max_ploidy = gt.nrow(); 

    for(int j=0; j<max_ploidy; j++) {
        for(int i=0; i<gens; i++){
            gt(j,_) = rcpp_mutate(gt(j,_), mu);
        }
    }
    return(gt);
}




