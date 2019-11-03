
#include "queue.h"

#define STRICT_R_HEADERS
#include <Rcpp.h>

/*
 * Logs statistics about the execution of the algorithm and dumps it to a file.
 * To turn off, pass verbosity <= 1
 */
NullLogger* logger;

//' R Interface to 'Certifiably Optimal RulE ListS (Corels)'
//'
//' Nothing here yet.
//' @title Corels interace
//' @return A constant bool for now
// [[Rcpp::export]]
bool corels() {
  return true;                  // more to fill in, naturally
}
