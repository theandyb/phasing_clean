# Rcpp
#include <Rcpp.h>
using namespace Rcpp;

//Return a vector of voted haplotypes
// Inputs are two vectors of haplotypes for each method(beagle, eagle, and shapeit)
// [[Rcpp::export]]
IntegerVector vote_function(IntegerVector beagle1, IntegerVector beagle2, 
    IntegerVector eagle1, IntegerVector eagle2, 
    IntegerVector shapeit1, IntegerVector shapeit2) {
  // Function implementation here
  List haplotypes;
  // Initialize empty vectors for the voted haplotypes
  IntegerVector voted_hap1(beagle1.size());
  //IntegerVector voted_hap2(beagle1.size());

  // Set current pointers to the first haplotype vectors
  IntegerVector* current_beagle = &beagle1;
  IntegerVector* current_eagle = &eagle1;
  IntegerVector* current_shapeit = &shapeit1;
  bool beagle_first = true;
  bool eagle_first = true;
  bool shapeit_first = true;


  // Loop through each position to perform voting
  for (int i = 0; i < beagle1.size(); i++) {
    // Check match between methods at current position for current pointer
    // First case: all match
    if ((*current_beagle)[i] == (*current_eagle)[i] && (*current_beagle)[i] == (*current_shapeit)[i]) {
      voted_hap1[i] = (*current_beagle)[i];
    } // second case: beagle and eagle match
    else if ((*current_beagle)[i] == (*current_eagle)[i]) {
      voted_hap1[i] = (*current_beagle)[i];
      // need to switch shapeit pointer
      if (shapeit_first) {
        current_shapeit = &shapeit2;
        shapeit_first = false;
      } else {
        current_shapeit = &shapeit1;
        shapeit_first = true;
      }
    } // third case: beagle and shapeit match
    else if ((*current_beagle)[i] == (*current_shapeit)[i]) {
      voted_hap1[i] = (*current_beagle)[i];
      // need to switch eagle pointer
      if (eagle_first) {
        current_eagle = &eagle2;
        eagle_first = false;
      } else {
        current_eagle = &eagle1;
        eagle_first = true;
      }
    } // fourth case: eagle and shapeit match
    else {
      voted_hap1[i] = (*current_eagle)[i];
      // need to switch beagle pointer
      if (beagle_first) {
        current_beagle = &beagle2;
        beagle_first = false;
      } else {
        current_beagle = &beagle1;
        beagle_first = true;
      }
    } 

  return voted_hap1;
}