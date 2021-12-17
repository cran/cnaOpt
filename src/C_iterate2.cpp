
#include <Rcpp.h>
#include "typedefs.h"
#include "C_iterate.h"

using namespace Rcpp;

typedef ListOf<NumericVector> dblList; 
typedef std::vector<NumericVector::iterator> dblIteratorList;
typedef std::pair<int, double> paired;


// // [[Rcpp::export]]
// IntegerVector C_order(NumericVector x, bool decreasing = false) {
//   NumericVector sorted = clone(x).sort(decreasing = decreasing);
// 	IntegerVector out=match(sorted, x);
// 	return out-1;
// }


// order() in Rcpp - code taken and adapted from here:
// https://stackoverflow.com/questions/21609934/ordering-permutation-in-rcpp-i-e-baseorder/26977061

bool cmp_second(const paired & left, const paired & right) {
    return left.second > right.second;
}

// [[Rcpp::export]]
Rcpp::IntegerVector C_order_desc(const Rcpp::NumericVector & x) {
    const size_t n = x.size();
    std::vector<paired> pairs;
    pairs.reserve(n);
    for(size_t i = 0; i < n; i++)
        pairs.push_back(std::make_pair(i, x(i)));
    std::sort(pairs.begin(), pairs.end(), cmp_second);
    Rcpp::IntegerVector result = Rcpp::no_init(n);
    for(size_t i = 0; i < n; i++)
        result(i) = pairs[i].first;
    return result;
}



/* R
library(microbenchmark)
x <- rnorm(1E5)
identical(order(-x), C_order_desc(x)+1L)
microbenchmark(
  order(-x),
  C_order_desc(x)
)
*/

// [[Rcpp::export]]
LogicalVector C_ccoKeep(const NumericVector con, const NumericVector cov){
	int n=con.size();
	LogicalVector out(n);
	if (n==0) return out;
	if (n==1){
		out(0)=true;
		return out;
	}
	IntegerVector ord=C_order_desc(con);
	int j=ord(0);
	double cummaxCov=cov(j);
	double lastCon=con(j);
	out(j)=true;
	for (int i=1; i<n; ++i){
		int j = ord(i);
		bool keep=(cov(j)>=cummaxCov);
		if ((cov(j)==cummaxCov) & (con(j)<lastCon)) keep=false;
		if (keep){
			out(j)=true;
			cummaxCov=cov(j);
			lastCon=con(j);
		}
		/*Rcpp::Rcout << "j=" << j << ", con=" << con(j) << ", cov=" << cov(j) << 
		 * ", keep=" << out(j) << ", cummaxCov=" << cummaxCov << std::endl; */
	}
	return out;
}

/*

n <- 50
x <- sample(as.numeric(0:10), n, T)/10
y <- sample(as.numeric(0:10), n, T)/10
keep <- C_ccoKeep(x, y) 
C_ccoKeep(x, y)
plot(x, y) 
points(x[keep], y[keep], col = "red", pch = 19)
  
*/


// [[Rcpp::export]]
List C_msubset(const dblList x, const LogicalVector s){
	int l=x.size();
	List out(l);
	for (int i=0; i<l; ++i){
		out[i]=x[i][s];
	}
	return out;
}

/*

x <- list(as.numeric(1:3), as.numeric(4:6), as.numeric(7:9))
s <- c(T, F, T)
C_msubset(x, s) 
C_msubset(C_msubset(x, s), c(F, T))

*/

// [[Rcpp::export]]
List C_getOptim(const dblList x){
	LogicalVector keep1=C_ccoKeep(x[0], x[1]);
	List x1=C_msubset(x, keep1);
	LogicalVector keep2=C_ccoKeep(x1[1], x1[0]);
	return C_msubset(x1, keep2);
}



// [[Rcpp::export]]
NumericVector C_append(const NumericVector x, const NumericVector y){
	int lx=x.size(), ly=y.size();
	NumericVector out(lx+ly);
	for (int i=0; i<lx; ++i){
		out(i)=x(i);
	}
	for (int i=0; i<ly; ++i){
		out(lx+i)=y(i);
	}	
	return(out);
}

// [[Rcpp::export]]
List C_mappend(const dblList x, const dblList y){
	int l=x.size();
	List out(l);
	for (int i=0; i<l; ++i){
		out[i]=C_append(x[i], y[i]);
	}
	return out;
}

/*

C_append(as.numeric(1:5), as.numeric(7:9))
C_mappend(
  list(as.numeric(1:3), as.numeric(6:4)), 
  list(as.numeric(4:6), as.numeric(3:1)))
C_mappend(
  list(as.numeric(1:3), as.numeric(6:4), pi), 
  list(as.numeric(4:6), as.numeric(3:1), exp(1))
 )
 
*/


// [[Rcpp::export]]
LogicalVector repTrueFalse(int len,int ntrue){
  LogicalVector out(len);
  for (int i=0; i<ntrue; ++i){
    out(i)=true;
  }
  return out;
}

// -----------------------------------------------------------------------------

// iteration using dblIteratorList
// [[Rcpp::export]]
List C_iterate2(dblList dx, dblList dminxy, double Sx_base, double Sy, 
                int blksize=10000, bool verbose = true){
  // initialize
  dblIteratorList dx_starts = getStarts(dx), dx_ends = getEnds(dx),
    dminxy_starts = getStarts(dminxy), dminxy_ends = getEnds(dminxy);
  dblIteratorList dx_it = dx_starts, dminxy_it = dminxy_starts;
  
  // Define output matrix
  NumericVector accum_con(blksize), accum_cov(blksize), accum_id(blksize);
  int pos = 0, changed = 0;
  long long int count = 0L;
  size_t n = dx.size();
  std::vector<double> cum_x(n), cum_minxy(n);
  dblList out=List::create(NumericVector(0), NumericVector(0), NumericVector(0));
  do {
    for (size_t ch = changed; ch < n; ++ch){
      if (ch == 0){
        cum_x.front() = *dx_it.front(); cum_minxy.front() = *dminxy_it.front();
      } else {
        cum_x[ch] = cum_x[ch-1] + *dx_it[ch]; 
        cum_minxy[ch] = cum_minxy[ch-1] + *dminxy_it[ch];
      }
    }
    double numerator = Sx_base + cum_minxy.back();
    double con = numerator / (Sx_base + cum_x.back());
    double cov = numerator / Sy;
    accum_con(pos) = con; accum_cov(pos) = cov; accum_id(pos) = static_cast<double>(count+1L);
    count++;
    pos++;
    if (pos==blksize){
      List accum=List::create(accum_con, accum_cov, accum_id);
      out=C_getOptim(C_mappend(out, accum));
      pos=0;
      if (verbose){
      	double n=static_cast<double>(ll(dx));
      	double pct=ceil(static_cast<double>(count+1L)/n*100);
    		Rcpp::Rcout << "ConCovOpt() progress: ~" << pct << "%\r" << std::flush;
    		}
    }
    increase(dx_it, changed, dx_starts, dx_ends);
    increase(dminxy_it, changed, dminxy_starts, dminxy_ends);
  } while(!finished(dx_it, dx_ends));
  if (pos>0){
    List accum=List::create(accum_con, accum_cov, accum_id);
    accum=C_msubset(accum, repTrueFalse(blksize, pos));
    out=C_getOptim(C_mappend(out, accum));
  }
  if (verbose) Rcpp::Rcout << "                           \r" << std::flush;
  return out;
}
