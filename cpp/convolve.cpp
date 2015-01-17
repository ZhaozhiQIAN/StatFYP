#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector convolveCpp(NumericVector a, NumericVector b) {
int na = a.size(), nb = b.size();
int nab = na + nb - 1;
NumericVector xab(nab);
for (int i = 0; i < na; i++)
for (int j = 0; j < nb; j++)
xab[i + j] += a[i] * b[j];
return xab;
}

// [[Rcpp::export]]
NumericVector ExtDiag(NumericMatrix a){
	int size = a.nrow();
	NumericVector vec(size);
	for (int i=0; i<size; ++i){
		vec[i] = a(i,i);
	}
	return vec;
}

// [[Rcpp::export]]
NumericVector CWeightGivenPi(NumericVector r1, NumericVector r2){

	int nobj = r2.size();
	int nrow = r1.size()/nobj;
	typedef NumericVector::iterator vec_iterator;
	vec_iterator itr1=r1.begin(),itr2=r2.begin();
	
	double L,I;
	NumericVector w(nrow*(nobj-1));
	vec_iterator itrw=w.begin();
	for (int k=0; k<nrow; ++k){
		for (int i=0; i<nobj; ++i){
			I=0;
			L=0;
			for (int j=0 ;j<nobj; ++j){
				if( ((itr1[i*nrow+k]>itr1[j*nrow+k]) && (itr2[i]<itr2[j])) || ((itr1[i*nrow+k]<itr1[j*nrow+k]) && (itr2[i]>itr2[j]))){
					++I;
				}
			}
			L= (itr1[i*nrow+k] + itr2[i] + I)/2;
			// to be optimized
			if(itr1[i*nrow+k] <= (L-1)){
				for (int p=itr1[i*nrow+k]-1; p<L-1;++p){
					itrw[p*nrow+k] = itrw[p*nrow+k]+0.5;
				}
			}
			if (itr2[i]<=(L-1)){
				for (int p=itr2[i]-1; p<L-1;++p){
					itrw[p*nrow+k] = itrw[p*nrow+k]+0.5;
				}
			}
		}
		// cumsum
		for (int i=1; i<nobj-1; ++i){
			itrw[i*nrow+k] = itrw[(i-1)*nrow+k] + itrw[i*nrow+k];
		}
	}
	return w;
}

// [[Rcpp::export]]
NumericVector CWeightGivenPiV2(NumericVector r1, NumericVector r2){

	int nobj = r2.size();
	int nrow = r1.size()/nobj;
	typedef NumericVector::iterator vec_iterator;
	vec_iterator itr1=r1.begin(),itr2=r2.begin();
	
	double I;
	double* pos1=new double[nobj];
	double* pos2=new double[nobj];
	double* L=new double[nobj];
	NumericVector w(nrow*(nobj-1));
	vec_iterator itrw=w.begin();
	for (int k=0; k<nrow; ++k){
		for (int i=0; i<nobj; ++i){
			I=0;
			for (int j=0 ;j<nobj; ++j){
				if( ((itr1[i*nrow+k]>itr1[j*nrow+k]) && (itr2[i]<itr2[j])) || ((itr1[i*nrow+k]<itr1[j*nrow+k]) && (itr2[i]>itr2[j]))){
					++I;
				}
			}
			pos1[i]=itr1[i*nrow+k];
			pos2[i]=itr2[i];
			L[i]= (pos1[i] + pos2[i] + I)/2;
		}
		for (int p=0; p<nobj; ++p){
			for (int q=0; q<nobj; ++q){
				if(pos1[q]-1 <= p && p< L[q]-1){
					itrw[p*nrow+k] += 0.5;
				}
				if(pos2[q]-1 <= p && p< L[q]-1){
					itrw[p*nrow+k] += 0.5;
				}
			}
		}
		// cumsum
		for (int i=1; i<nobj-1; ++i){
			itrw[i*nrow+k] = itrw[(i-1)*nrow+k] + itrw[i*nrow+k];
		}
	}
	return w;
}












