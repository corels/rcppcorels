
#include "queue.h"

#include <Rcpp.h>

/*
 * Logs statistics about the execution of the algorithm and dumps it to a file.
 * To turn off, pass verbosity <= 1
 */
NullLogger* logger;

// [[Rcpp::export]]
bool corels() {
  return true;                  // more to fill in, naturally
}
