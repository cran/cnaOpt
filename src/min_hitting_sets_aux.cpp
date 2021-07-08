
#include <Rcpp.h>
#include "typedefs.h"

using namespace Rcpp;

bool C_intersects_with_all(CharacterVector x, charList y){
	bool out=true;
	int niter=y.size();
	for (int i=0; i<niter; ++i){
		CharacterVector intersection=intersect(x, y[i]);
		if (intersection.size() == 0){
			out=false;
			break;
		}
	}
	return(out);
}

// [[Rcpp::export]]
LogicalVector C_m_intersects_with_all(charList x, charList y){
	int n=x.size();
	LogicalVector out(n);
	for (int i=0; i<n; ++i){
		out[i] = C_intersects_with_all(x[i], y);
	}
	return(out);
}


bool C_contains_one_of(CharacterVector x, charList y){
	bool out=false;
	int niter=y.size();
	for (int i=0; i<niter; ++i){
		CharacterVector setdifference = setdiff(y[i], x);
		if (setdifference.size() == 0){
			out=true;
			break;
		}
	}
	return(out);
}


// [[Rcpp::export]]
LogicalVector C_m_contains_one_of(charList x, charList y){
	int n=x.size();
	LogicalVector out(n);
	for (int i=0; i<n; ++i){
		out[i] = C_contains_one_of(x[i], y);
	}
	return(out);
}

