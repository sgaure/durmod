#include <cstdlib>
#include <utility>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <Rcpp.h>
#include <omp.h>
#include <string.h>
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;

enum Timing {exact, interval, logit};

typedef struct {
  int *val;
  int nlevels;
  double *x;
} FACTOR;

typedef struct {
  double *par;
  int nlevels;
} FACTORPAR;


static bool debug = false;
// [[Rcpp::export]]
void cdebug(bool dodebug) {
  debug=dodebug;
}

// return log(sum(exp(x+y)))
inline double logsumofexp(const int n,const double *x, const double *y) {
  if(n < 1) stop("logsumofexp called with n < 1 (%d)",n);
  double cur = x[0]+y[0];
  for(int i = 1; i < n; i++) cur = R::logspace_add(cur,x[i]+y[i]);
  return cur;
}

// log(1+sum(x))
inline double log1sumx(int n, const double *x) {
  double lsum = 0;
  for(int i = 0; i < n; i++) lsum = R::logspace_add(lsum,x[i]);
  return lsum;
}



// convert parametrized probabilities Pi = exp(ai)/(1+sum(exp(ai))) to log probabilities
// a2logp <- function(a) {b <- c(0,a); logp <- b - logsumofexp(b,0); ifelse(is.na(logp),0,logp)}
void a2logp(const int n, const double *a, double *logp) {
  if(n < 0) stop("a2logp called with negative n (%d)",n);
  
  double lsum = log1sumx(n, a);
  logp[0] = -lsum;
  for(int i = 0; i < n; i++) logp[i+1] = a[i]-lsum;
  return;
}


// Taken from Martin Maechler https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
// "Accurately Computing log(1-exp(-a)), Assessed by the Rmpfr package"
// compute log(1-exp(-x))
inline double log1mexp(double x) {
  return (x <= M_LN2) ? log(-expm1(-x)) : log1p(-exp(-x));
}

/*
// log(1+exp(x))
inline double log1pexp(double x) {
  if(x <= -37) return exp(x);
  if(x <= 18) return log1p(exp(x));
  if(x <= 33.3) return x+exp(-x);
  return x;
}
*/
// compute exp(-x)/(1-exp(-x))
// by logarithms exp(-x - log(1-exp(-x)))
inline double expfrac(double x,double extra = 0.0) {
  return exp(-x - log1mexp(x) + extra);
}


//      computeloglik(d[i], lh, mup, npoints, llspell)
inline void obsloglik(const int tr, const Timing timing, const double *lh, double dur,
		      double **mup, const int npoints, int transitions,
		      const bool *riskmask,
		      double *llspell) {

  // If there is a transition, add the loghazard
  // divide by the survival prob up until now, i.e. subtract its log
  for(int j = 0; j < npoints; j++) {
    double logsumhaz = -DBL_MAX;
    for(int t = 0; t < transitions; t++) {
      if(riskmask && !riskmask[t]) continue;
      logsumhaz = R::logspace_add(logsumhaz, lh[t]+mup[t][j]);
    }
    switch(timing) {
    case exact:
    case interval:
      llspell[j] -= dur*exp(logsumhaz);
      break;
    case logit:
      llspell[j] -= logsumhaz;
      break;
    }
    if(tr > 0) {
      switch(timing) {
      case interval:
	llspell[j] += log1mexp(dur*exp(logsumhaz)) - logsumhaz;
      case exact:
      case logit:
	llspell[j] += lh[tr-1] + mup[tr-1][j];
	break;
      }
    }
  }
}

// compute gradient
inline void gobsloglik(const int tr, const Timing timing, const double *lh, double dur, int obs,
		       double **mup, const int npoints, int transitions,
		       int npars, int *nfacs,
		       const bool *riskmask,
		       int *nvars,
		       double **matp,
		       FACTOR **factors,
		       double *dllspell) {
  int inipos = 0;
  const int t = tr-1;
  for(int s = 0; s < t; s++) inipos += nvars[s] + nfacs[s] + npoints;
  const double *tmat = (t>=0) ? &matp[t][obs*nvars[t]] : 0;

  for(int j = 0; j < npoints; j++) {
    double sumhaz = 0;
    double logsumhaz = -DBL_MAX;
    if(timing == logit || (timing == interval && t >= 0)) {
      for(int tt = 0; tt < transitions; tt++) {
	if(riskmask && !riskmask[tt]) continue;
	logsumhaz = R::logspace_add(logsumhaz, lh[tt] + mup[tt][j]);
      }
      sumhaz = exp(logsumhaz);
    }

    double *dll = &dllspell[j*npars];
    int pos = 0;
    for(int tt = 0; tt < transitions; tt++) {
      if(riskmask && !riskmask[tt]) {pos += nvars[tt]+nfacs[tt]+npoints;continue;}
      const double haz = exp(lh[tt] + mup[tt][j]);
      const double *mat = &matp[tt][obs*nvars[tt]];

      for(int k = 0; k < nvars[tt]; k++) {
	switch(timing) {
	case exact:
	case interval:
	  dll[pos++] -= dur*haz*mat[k];
	  break;
	case logit:
	  dll[pos++] -= haz*mat[k]/sumhaz;
	  break;
	}
      }
      const FACTOR *fac = factors[tt];
      // loop through the factors for this transition
      for(int k = 0; k < nfacs[tt]; k++) {
	const int fval = fac[k].val[obs];
	if(fval <= 0) {pos += fac[k].nlevels; continue;}
	const double f = (fac[k].x != 0) ? fac[k].x[obs] : 1.0;
	switch(timing) {
	case exact:
	case interval:
	  dll[pos + fval-1] -= dur*haz*f;
	  break;
	case logit:
	  dll[pos + fval-1] -= haz*f/sumhaz;
	  break;
	}
	pos += fac[k].nlevels;
      }
      // The mu, the covariate is 1
     
      switch(timing) {
      case exact:
      case interval:
	dll[pos+j] -= dur*haz;
	break;
      case logit:
	dll[pos+j] -= haz/sumhaz;
	break;
      }
      pos += npoints;
    }
    if(t >= 0) {
      //      double *dll = &dllspell[j*npars];
      int pos = inipos;
      // complicated expression for interval timing, the derivative of llspell wrt to a parameter
      // With H the sum of hazards, L the log hazard of transition t, and d duration, we compute
      // the derivative w.r.t x of the log likelihood
      // log(1-exp(-d H)) - log(H) + L
      
      double funnyexpr;
      if(timing == interval)
	funnyexpr = dur*exp(lh[t]+mup[t][j] - dur*sumhaz)/expm1(-dur*sumhaz) +
	  expm1(lh[t]+mup[t][j]-logsumhaz);

      for(int k = 0; k < nvars[t]; k++) {
	switch(timing) {
	case exact:
	case logit:
	  dll[pos++] += tmat[k];
	  break;
	case interval:
	  dll[pos++] -= tmat[k]*funnyexpr;
	  break;
	}
      }
      
      const FACTOR *fac = factors[t];
      for(int jj = 0; jj < nfacs[t]; jj++) {
	const int fval = fac[jj].val[obs];
	if(fval <= 0) {pos += fac[jj].nlevels; continue;};  // skip NA-levels, i.e. reference
	double *x = fac[jj].x;
	double f = (x != 0) ? x[obs] : 1.0;
	switch(timing) {
	case exact:
	case logit:
	  dll[pos + fval-1] += f;
	  break;
	case interval:
	  dll[pos + fval-1] -= f*funnyexpr;
	  break;
	}
	pos += fac[jj].nlevels;
      }
      
      // and for the mus
      switch(timing) {
      case exact:
      case logit:
	dll[pos+j] += 1;
	break;
      case interval:
	dll[pos+j] -= funnyexpr;
	break;
      }
    }
  }
}

inline void updategradient(int npoints, double *dllspell, double *llspell, double *logprobs, double ll,
		     int transitions, int npars, int *nvars, int *faclevels, const double *pargs, 
		     int totalpars, double *spellgrad, double *grad) {

  (void) memset(spellgrad, 0, totalpars*sizeof(double));

  // compute gradient
  // we have the gradient of each llspell component in dllspell
  // Now, ll = log(sum(p_j exp(lh_j)))
  // where j ranges over 1..masspoints
  // We have lh_j in llspell
  // we have d lh_j / dx in dllspell
  // we have d ll / dx = 1/sum(p_j exp(lh_j)) * d sum(p_j exp(lh_j)) / dx.
  //
  // For x not a probability parameter, i.e. only occuring in lh_j
  // we have d sum(p_j exp(lh_j)) / dx = sum(p_j exp(lh_j) dlh_j/dx).
  // 
  // For the probability parameters y (which do not occur in each of the lh_j) we have
  // d sum(p_j exp(lh_j)) / dy = sum(dp_j/dy exp(lh_j))
  
  // then compute d sum(p_j exp(lh_j)) / dx.
  // First for the ordinary covariates where it equals sum(p_j exp(lh_j) dlh_j/dx)

  for(int j = 0; j < npoints; j++) {
    const double *dll = &dllspell[j*npars];
    const double scale = exp(logprobs[j] + llspell[j] - ll);
    
    int pos = 0;
    for(int t = 0; t < transitions; t++) {
      for(int k = 0; k < nvars[t]+faclevels[t]; k++) {
	spellgrad[pos] += scale * dll[pos];
	pos++;
      }
      // the mu for this masspoint in this transition
      spellgrad[pos+j] += scale*dll[pos+j];
      pos += npoints; // next transition
    }
    
    // The probability-parameters ak occur in all the probabilities
    // compute dPj / dak
    // The a's are in pargs
    // let ld = log(1+sum(exp(pargs))), can be moved out
    // Stick to logs, 
    double ld = log1sumx(npoints-1, pargs);
    const double lscale = llspell[j] - ll;
    for(int k = 0; k < npoints-1; k++) {
      const double ak = pargs[k];
      double dPdak;

      if(j == 0) {
	// for j=0, special case, it's dP0/dak
	dPdak = -exp(ak-2*ld + lscale);
      } else if(j == k+1) {
	// dPk/dak
	dPdak = exp(ak-ld + lscale) - exp(2*(ak-ld) + lscale);
      } else {
	// dPj/dak
	dPdak = -exp(logprobs[k+1] + logprobs[j] + lscale);
      }

      spellgrad[npars+k] += dPdak; 
    }
  }
  // update global gradient with spell gradient
  for(int k = 0; k < totalpars; k++) grad[k] += spellgrad[k];

}

inline void  updatefisher(int *gradfill, int fishblock, int totalpars, double *gradblock, 
			  int *nonzero, double *spellgrad, double *fisher) {
  // When computing the fisher matrix, we have individual spell
  // gradients in the gradblock matrix, we should add this
  // spell's gradient as a column.  When the gradblock matrix is
  // full, we dsyrk it into the fisher matrix. We do this in a
  // critical region, since all the threads use the same
  // gradblock and fisher matrix. This has the additional
  // bonus that only one thread runs the dsyrk at any time, so
  // we can use a parallel blas. The reason we collect gradients
  // in gradblock is that dsyrk is a lot faster than dsyr (rank
  // 1 update).


  if(*gradfill == fishblock) {
    // dsyrk the gradblock into fisher
    //       subroutine dsyrk (UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC)
    // C = alpha A * A' + beta C
    
    // how many rows are nonzero?
    int nnrank = 0;
    for(int j = 0; j < totalpars; j++) {
      bool nonz = false;
      for(int k = 0; k < fishblock && !nonz; k++) {
	nonz |= (gradblock[k*totalpars + j] != 0.0);
      }
      if(nonz) nonzero[nnrank++] = j;
    }
    if(4*nnrank < 3*totalpars) {
      // Less than some limit, we dsyrk a smaller one
      
      double *smallblock = new double[fishblock*nnrank];
      for(int b = 0; b < fishblock; b++) {
	for(int k = 0; k < nnrank; k++) {
	  smallblock[b*nnrank + k] = gradblock[b*totalpars + nonzero[k]];
	}
      }
      // dsyrk it into smallfish
      const double alpha=-1, beta=0;
      
      double *smallfish = new double[nnrank*nnrank];
      F77_CALL(dsyrk)("U","N",&nnrank,&fishblock,&alpha,smallblock,&nnrank,&beta,
		      smallfish, &nnrank);
      delete [] smallblock;
      // update fisher from smallfish
      for(int j = 0; j < nnrank; j++) {
	for(int k = 0; k <= j; k++) {
	  fisher[nonzero[j]*totalpars + nonzero[k]] += smallfish[j*nnrank + k];
	}
      }
      delete [] smallfish;
    } else {
      const double alpha = -1, beta = 1;
      F77_CALL(dsyrk)("U","N", &totalpars, &fishblock, &alpha, gradblock, 
		      &totalpars, &beta, fisher, &totalpars);
    }
    *gradfill = 0;
  }
  // append the spell gradient to the block
  for(int k = 0; k < totalpars; k++) gradblock[*gradfill * totalpars + k] = spellgrad[k];
  (*gradfill)++;
}

// [[Rcpp::export]]
NumericVector cloglik(List spec, List pset, 
		      const bool gdiff=false, const bool dogradient=false, const bool dofisher=false,
		      const int nthreads=1) {
  const IntegerVector d = as<IntegerVector>(spec.attr("d"));
  const IntegerVector id = as<IntegerVector>(spec.attr("id"));
  const NumericVector duration = as<NumericVector>(spec.attr("duration"));
  const NumericVector spellidx = as<NumericVector>(spec.attr("spellidx"));
  const CharacterVector ctiming = as<CharacterVector>(spec.attr("timing"));
  const Timing timing = (ctiming[0] == "exact") ? exact : (ctiming[0] == "interval" ? interval : logit);

  const int nspells = spellidx.size() - 1;
  List parset = as<List>(pset["parset"]);
  const double *pargs = REAL(as<NumericVector>(pset["pargs"]));
  const int npoints = 1 + as<NumericVector>(pset["pargs"]).size();
  double *logprobs = new double[npoints];

  a2logp(npoints-1,pargs,logprobs);
  const int N = d.size();
  const int transitions = parset.size();

  List risklist = as<List>(spec.attr("riskset"));
  const int nrisks = risklist.size();
  IntegerVector state = as<IntegerVector>(spec.attr("state"));

  bool *riskmasks = new bool[nrisks*transitions]();
  for(int i = 0; i < nrisks; i++) {
    IntegerVector v = as<IntegerVector>(risklist[i]);
    for(int t = 0; t < v.size(); t++) {
      riskmasks[i*transitions + v[t]-1] = true;
    }
  }
  
  // pointers into the data

  double **matp = new double*[transitions];

  int *nvars = new int[transitions];

  // number of factors in each transition

  int *nfacs = new int[transitions];
  // pointers to factor lists
  FACTOR **factors = new FACTOR*[transitions];
  for(int i = 0; i < transitions; i++) {
    NumericMatrix smat = as<NumericMatrix>(as<List>(spec[i])["mat"]);
    matp[i] = REAL(smat);
    nvars[i] = smat.nrow();

    List facs = as<List>(spec[i])["faclist"];
    nfacs[i] = facs.size();
    //    printf("nfacs[%d] = %ld\n",i,nfacs[i]);
    if(nfacs[i] == 0) continue;
    factors[i] = new FACTOR[nfacs[i]];
    for(int j = 0; j < nfacs[i]; j++) {
      IntegerVector fac = as<IntegerVector>(facs[j]);
      factors[i][j].val = INTEGER(fac);
      factors[i][j].nlevels = as<CharacterVector>(fac.attr("levels")).size();
      NumericVector xv = as<NumericVector>(fac.attr("x"));
      if(xv.size() > 0) {
	if(xv.size() != N) stop("Factor interaction term length(%d) does not match dataset(%d)",
				xv.size(), N);
	factors[i][j].x = REAL(xv);
      } else {
	factors[i][j].x = 0;
      }
    }
  }

  // pointers into the parameters

  double **betap = new double*[transitions];
  double **mup = new double*[transitions];


  FACTORPAR **facpars = new FACTORPAR*[transitions];
  int *faclevels = new int[transitions];
  for(int i = 0; i < transitions; i++) {
    betap[i] = REAL(as<List>(parset[i])["pars"]);
    mup[i] = REAL(as<List>(parset[i])["mu"]);
    List Rfacs = as<List>(parset[i])["facs"];
    if(nfacs[i] != Rfacs.size()) stop("number of factor parameters(%d) does no match data(%d)",
				      Rfacs.size(),nfacs[i]);
    faclevels[i] = 0;
    if(nfacs[i] > 0) {
      facpars[i] = new FACTORPAR[Rfacs.size()];
      for(int j = 0; j < nfacs[i]; j++) {
	facpars[i][j].par = REAL(as<NumericVector>(Rfacs[j]));
	facpars[i][j].nlevels = as<NumericVector>(Rfacs[j]).size();
	faclevels[i] += facpars[i][j].nlevels;
      }
    }
  }

  const int fishblock = 128 ;
  const bool dograd = dogradient ? true : dofisher;

  // find the number of parameters
  int npars = 0;
  for(int t = 0; t < transitions; t++) npars += nvars[t] + faclevels[t];
// and mu's
  npars += npoints*transitions;

  const int totalpars = npars+npoints-1;

  double LL = 0.0;
  const int gradsize = dograd ? totalpars : 1;
  double *grad = new double[gradsize]();

  double *gradblock = 0;
  if(dofisher) gradblock = new double[fishblock*gradsize];
  int gradfill = 0;

  // Our fisher matrix is global, and updated inside a critial region
  // We could do a reduction(+) on it instead, but that may consume too much memory
  NumericMatrix *retfisher;
  double *fisher;
  if(dofisher) {
    retfisher = new NumericMatrix(gradsize,gradsize);
    fisher = REAL(*retfisher);
  }

  // Then some thread private storage which are static in the parallel for below.
  // We don't want to allocate in the loop.
  static double *llspell, *lh, *dllspell, *spellgrad;  // must be static to allocate thread private
  static int *nonzero;
#pragma omp threadprivate(llspell,dllspell, lh, spellgrad, nonzero)
  // Allocate it
#pragma omp parallel num_threads(nthreads)
  {  
    llspell = new double[npoints];
    lh = new double[transitions];
    if(dograd) {
      dllspell = new double[npoints*npars];
      spellgrad = new double[totalpars];
    }
    if(dofisher) nonzero = new int[totalpars];

  }

  // Remember not to use any R-functions (or allocate Rcpp storage) inside the parallel region.

#pragma omp parallel for reduction(+:grad[:gradsize], LL) num_threads(nthreads) schedule(guided)
  for(int spellno = 0; spellno < nspells; spellno++) {
    memset(llspell, 0, npoints*sizeof(double));
    if(dograd) memset(dllspell, 0, npoints*npars*sizeof(double));
    for(int i = spellidx[spellno]; i < spellidx[spellno+1]; i++) {
      
      // compute the log hazard for this observation for each masspoint and transition
      // loop through the mass points. For each find the hazard sum
      // fill in the loghazards in the lh-array
      const bool *riskmask = nrisks>0 ? &riskmasks[(state[i]-1)*transitions] : 0;
      for(int t = 0; t < transitions; t++) {
	lh[t] = 0.0;
	if(riskmask && !riskmask[t]) continue;  // should this be here?
	const double *mat = &matp[t][i*nvars[t]];
	const double *beta = betap[t];
	for(int k = 0; k < nvars[t]; k++) {
	  lh[t] += mat[k]*beta[k];
	}
	const FACTOR *fac = factors[t];
	const FACTORPAR *fpar = facpars[t];
	for(int j = 0; j < nfacs[t]; j++) {
	  const int fval = fac[j].val[i];
	  if(fval <= 0) continue;  // skip NA-levels, i.e. reference
	  double *x = fac[j].x;
	  if(x != 0) lh[t] += fpar[j].par[fval-1] * x[i]; else lh[t] += fpar[j].par[fval-1];
	}
      }
      

      // update llspell with the observation log likelihood 
      obsloglik(d[i], timing, lh, duration[i], mup, npoints, transitions, riskmask, llspell);
      // update dllspell with the gradient of the observation log likelihood
      if(dograd) gobsloglik(d[i], timing, lh, duration[i], i, mup, npoints, transitions, npars, nfacs,
			    riskmask, nvars, matp, factors, dllspell);
    }
    
    // We have collected the loglikelihood of a spell, one for each masspoint
    // integrate it with the probabilities
    double ll;
    if(gdiff) {
      ll = logsumofexp(npoints-1,llspell,logprobs);
      LL += exp(llspell[npoints-1] - ll);
    } else {
      // compute the log likelihood
      ll = logsumofexp(npoints,llspell,logprobs);
      LL += ll;
      
      if(dograd) {
	// compute the spell gradient, update the gradient
	updategradient(npoints, dllspell, llspell, logprobs, ll,
		 transitions, npars, nvars, faclevels, pargs, totalpars, spellgrad, grad);
	if(dofisher) {
	  // update the fisher matrix from the spellgrad
	  // we use a global fisher matrix, no omp reduction, so do it in a critical section
#pragma omp critical
	  updatefisher(&gradfill, fishblock, totalpars, gradblock, nonzero, spellgrad, fisher);
	}
      }
    } 
  }  // end of spell loop

  // Deallocate thread private storage
#pragma omp parallel num_threads(nthreads)
  {  
    delete [] llspell;
    delete [] lh;
    if(dograd) {
      delete [] dllspell;
      delete [] spellgrad;
    }
    if(dofisher) delete [] nonzero;
  }


  // if anything remains in the gradient blocks, dsyrk it into the fisher matrix
  if(dofisher && gradfill > 0) {
    const double alpha = -1, beta = 1;
    const int N = totalpars, K=gradfill;
    F77_CALL(dsyrk)("U","N", &N, &K, &alpha, gradblock, &N, &beta, fisher, &N);
 
  }

  // deallocate other storage
  if(dofisher) delete [] gradblock;

  delete [] logprobs;
  delete [] riskmasks;
  delete [] matp;
  delete [] nvars;
  for(int t = 0; t < transitions; t++) {
    if(nfacs[t] > 0) delete [] factors[t];
    if(nfacs[t] > 0) delete [] facpars[t];
  }
  delete [] factors;
  delete [] nfacs;
  delete [] betap;
  delete [] mup;
  delete [] facpars;
  delete [] faclevels;

  // Then set up the return value
  NumericVector ret = NumericVector::create(LL);
  if(dograd) {
    NumericVector retgrad(totalpars);
    (void) memcpy(REAL(retgrad), grad, gradsize*sizeof(double));
    ret.attr("gradient") = retgrad;
    delete [] grad;
  }
  if(dofisher)  {
    // fill the lower half of the fisher matrix
    // read consecutively, write with stride, that's typically faster
    for(int i = 0; i < totalpars; i++) {
      for(int j = 0; j < i; j++) {
	fisher[j*totalpars + i] = fisher[i*totalpars + j];
      }
    }
    ret.attr("fisher") = *retfisher;
  }
  return ret;
}
