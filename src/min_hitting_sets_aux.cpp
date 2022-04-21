
#include <Rcpp.h>
#include <R_ext/Utils.h>
#include "typedefs.h"

using namespace Rcpp;

IntegerVector combn_initialize_indices(const int m){
	IntegerVector ind(m);
	for (int i=0; i<m; ++i){
		ind(i)=i;
	}
	return(ind);
}

void combn_increase(IntegerVector ind, const int n){
	int m=ind.size();
	int pos=m-1;
	for (int i=m-1; i>=0; --i){
		if (ind(i)<n-m+i){
			break;
		} else {
			pos=pos-1;
		}
	}
	if (pos<0){
		ind(0)=-1;
		return;
	}
	ind(pos)=ind(pos)+1;
	if (pos<m-1){
		for (int i=pos+1; i<m; ++i){
			ind(i)=ind(i-1)+1;
		}
	}
}

// [[Rcpp::export]]
List resize(const List x, const int n){
	if (n==0) return(List(0));
  int oldsize = x.size();
  List y(n);
  if (oldsize>n) oldsize=n;
  for( int i=0; i<oldsize; i++){
  	y[i] = x[i];
  }
  return y;
}

bool C_intersectsWithAll(const IntegerVector rows, const LogicalMatrix m) {
	int nc=m.ncol(), nr=rows.size();
	bool out=true;
  for (int j=0; j<nc; ++j){
  	bool found=false;
  	for (int i=0; i<nr; ++i){
  		int r=rows(i);
  		if (m(r, j)){
  			found=true;
  			break;
  		}
  	}
		if (!found){
			out=false;
			break;
	  }
  }
  return(out);
}

bool C_containsSolution(const IntegerVector rows, 
                        const LogicalMatrix m, 
                        const IntegerVector solLengths) {
	int nsol=m.ncol(), nr=rows.size();
	bool out=false;
  for (int s=0; s<nsol; ++s){
  	int count=0;
  	for (int i=0; i<nr; ++i){
  		int r=rows(i);
  		if (m(r, s)){
  			count++;
  		} else {
  		}
  	}
  	if (count == solLengths(s)){
  		out=true;
  		break;
  	}
  }
  return(out);
}


// [[Rcpp::export]]
List C_mhs_iteration(const int m, const LogicalMatrix m_x, const LogicalMatrix m_sol) {
	int n=m_x.nrow(), outpos=0;
	IntegerVector solLengths(0);
	if (m_sol.ncol()>0) solLengths=colSums(m_sol);
  List out(100);
  IntegerVector ind=combn_initialize_indices(m);
  while(ind(0)>=0){
  	if ((solLengths.size()>0) && C_containsSolution(ind, m_sol, solLengths)){
  		combn_increase(ind, n);
  		continue;
  	}
  	if (C_intersectsWithAll(ind, m_x)){
  		if (outpos>=out.size()){
  			out=resize(out, outpos+100);
  		}
  		out[outpos]=clone(ind);
  		++outpos;
  	}
  	combn_increase(ind, n);
  	R_CheckUserInterrupt();
  }
  out=resize(out, outpos);
  return(out);
}

