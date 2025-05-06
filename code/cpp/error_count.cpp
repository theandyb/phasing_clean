// [[Rcpp::depends(RcppArmadillo)]]

#include "RcppArmadillo.h"
#include <vector>

// [[Rcpp::export]]
Rcpp::List switch_flip_counter(const arma::vec& true1,
                               const arma::vec& true2,
                               const arma::vec& test1,
                               const arma::vec& test2){
  size_t n = true1.n_elem;
  if(true2.n_elem != n || test1.n_elem != n || test2.n_elem != n){
    Rcpp::stop("Input vectors must have the same length.");
  }
  bool swap = false;
  int switch_count, flip_count;
  flip_count = 0;
  switch_count = 0;
  bool prev_switch = true;

  for(int i = 0; i < n; ++i){
    if(swap){
      // do we have a match on the current haplotype?
      if(test2(i) == true1(i)){
        prev_switch = false;
      }
      //A mismatch is either due to a switch or flip
      else{
        if(prev_switch){
          // we have two switches in a row, which means a flip!
          flip_count++;
          switch_count--;
          prev_switch = false;
        } else{
          // we might have a switch
          switch_count++;
          prev_switch = true;
        }
        // move back to the other haplotype for the next position
        swap = false;
      }
    }
    else {
      if(test1(i) == true1(i)){
        prev_switch = false;
      }
      else {
        if(prev_switch){
          flip_count++;
          switch_count--;
          prev_switch = false;
        }
        else{
          switch_count++;
          prev_switch = true;
        }
        swap = true;
      }
    }
  }

  return Rcpp::List::create(Rcpp::Named("switches") = switch_count,
                            Rcpp::Named("flips") = flip_count);
}

// [[Rcpp::export]]
Rcpp::List switch_flip(const arma::ivec& end_pos,
                       const arma::ivec& start_pos){
  std::vector<int> flip_pos;
  std::vector<int> switch_pos;

  bool in_progress = false;
  int t = end_pos.n_elem;

  for(int i = 0; i < (t - 1); ++i){
    if(end_pos[i] == start_pos[i+1]){
      if(!in_progress){
        in_progress = true;
        flip_pos.push_back(end_pos[i]);
      } else{
        // position is part of latter half of a flip followed by either a switch or a flip
        in_progress = false;
      }
    }
    else if(in_progress){
      // second half of a flip
      in_progress = false;
    }
    else{
      switch_pos.push_back(end_pos[i]);
    }
  }

  return Rcpp::List::create(Rcpp::Named("switches") = switch_pos,
                            Rcpp::Named("flips") = flip_pos);
}
