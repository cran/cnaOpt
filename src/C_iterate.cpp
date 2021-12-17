
#include <Rcpp.h>
#include "typedefs.h"

using namespace Rcpp;

// show() for dblIteratorList & std::vector<double>
// void show(dblIteratorList x){
//   for (size_t i=0; i<x.size(); ++i){
//     NumericVector::iterator xi = x[i];
//     Rcout << *xi << " ";
//   }
//   Rcout << std::endl;
// }
// void show(std::vector<double> x){
//   for (size_t i=0; i<x.size(); ++i){
//     Rcout << x[i] << " ";
//   }
//   Rcout << std::endl;
// }

// getStarts and getEnds for dblIteratorList
dblIteratorList getStarts(dblList x){
  dblIteratorList out;
  int l = x.size();
  for (int i=0; i<l; ++i){
    out.push_back(x[i].begin());
  }
  return out;
}
dblIteratorList getEnds(dblList x){
  dblIteratorList out;
  int l = x.size();
  for (int i=0; i<l; ++i){
    out.push_back(x[i].end());
  }
  return out;
}

// increase for dblIteratorList
void increase(dblIteratorList& ii, int& changed,
              const dblIteratorList& _starts_, 
              const dblIteratorList& _ends_, 
              int pos = -1){
  if (pos == -1) pos = ii.size() - 1;
  ++ii[pos];
  changed = pos;
  if (ii[pos] == _ends_[pos]){
    if (pos > 0){
      ii[pos] = _starts_[pos];
    increase(ii, changed, _starts_, _ends_, pos = pos-1);
    }
  }
}

bool finished(dblIteratorList ii, 
              const dblIteratorList _ends_){
  return (ii.front() == _ends_.front());
};

// Test function
//// [[Rcpp::export]]
// void test(dblList x){
//   dblIteratorList _starts_ = getStarts(x);
//   dblIteratorList _ends_ = getEnds(x);
//   dblIteratorList ii = _starts_;
//   int changed = 0;
//   do {
//     for (int ch=changed; ch < x.size(); ++ch){
//       Rcout << "      > Inserting "<< *ii[ch] <<
//         " at position " << ch << std::endl;
//     }
//     show(ii);
//     // more actions here if required...
//     increase(ii, changed, _starts_, _ends_);
//   } while(!finished(ii, _ends_));
// }


/* R
test(list(1:2, 5:9))
test(list(4:2 * 2, 4:8, c(55, 44, 66)))
test(list(1, 2, 3:4, 5, 6:7))
*/

// double sum(dblIteratorList x){
//   double out = 0;
//   for (unsigned int i=0; i<x.size(); ++i){
//     out += *x[i];
//   }
//   return out;
// }

// Aux fun dblList -> easier using std::accumulate!
long long int ll(dblList x){
  long long int out = 1;
  for (int i=0; i<x.size(); ++i){
    out *= static_cast<long long int>(x[i].size());
  }
  return out;
}

