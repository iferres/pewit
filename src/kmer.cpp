#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
CharacterVector compute_kmers( const std::string seq, const int k ) {

  const int len = seq.length();

  const int n1 = k - 1;

  CharacterVector result( len - n1 );

  int j = 1;

  std::string tmp = seq.substr( 0, k );

  result[0] = tmp;

  for ( int i = k; i < len; i++ ){

    tmp.erase( 0, 1 );

    tmp.push_back( seq[i] );

    result[j] = tmp;

    j++;
  }

  return result;
}
