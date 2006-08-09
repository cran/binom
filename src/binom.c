#include <math.h>

#define DBL_EPSILON 2.220446e-016
#define DBL_MIN 1e-14

double lgamma(double c) {
  long j;
  double x,y,tmp;
  double ser = 1.000000000190015;
  const double cof[6] = {
    76.18009172947146,
    -86.50532032941677,
    24.01409824083091,
    -1.231739572450155 ,
    0.1208650973866179e-2,
    -0.5395239384953e-5};
  y = x = c;
  tmp = x + 5.5 - (x + 0.5) * log(x + 5.5);
  for(j=0;j<=5;j++)
    ser += (cof[j] / ++y);
  return log(2.5066282746310005 * ser / x) - tmp;
}

inline double lfactorial(long n)
{ return lgamma(n+1); }

inline double expm1(double x) /* exp(x)-1 */
{ return exp(x)-1.0; }

inline double log1p(double x) /* log(1+x) */
{ return fabs(x)<DBL_EPSILON?(x-0.5*x*x):log(1+x); }

inline double lbeta(double a,double b) /* log of the beta function  */
{ return lgamma(a)+lgamma(b)-lgamma(a+b); }

inline double dbeta(double x,double a,double b)
{ return exp(-lbeta(a,b)+(a-1)*log(x)+(b-1)*log(1-x)); }

double pbeta(double x,double pin,double qin,long upper) {
  // double pbeta_raw(double x,double pin,double qin,bool upper);
  // ALGORITHM AS 63 APPL. STATIST. VOL.32, NO.1
  // Bosten and Battiste (1974).
  // Remark on Algorithm 179, CACM 17, p153, (1974).
  // Computes P(Beta>x)
  double ans, c, finsum, p, ps, p1, q, term, xb, xi, y;
  int n, i, ib;
  long swap_tail;  
  const double eps = 0.5*DBL_EPSILON;
  const double sml = DBL_MIN;
  const double lneps = log(eps);
  const double lnsml = log(sml);
  /* swap tails if x is greater than the mean */
  if(pin / (pin + qin) < x) {
    swap_tail = 0;
    y = 1 - x;
    p = qin;
    q = pin;
  } else {
    swap_tail = 1;
    y = x;
    p = pin;
    q = qin;
  }
  if((p + q) * y / (p + 1) < eps) {
    /* tail approximation */
    xb = p * log(fmax(y, sml)) - log(p) - lbeta(p, q);
    if(xb > lnsml && y != 0)
      ans = (swap_tail == upper) ? exp(xb) : -expm1(xb);
    else
      ans = (swap_tail == upper) ? 0. : 1.;
  } else {
    /* evaluate the infinite sum first.  term will equal */
    /* y^p / beta(ps, p) * (1 - ps)-sub-i * y^i / fac(i) */
    ps = q - floor(q);
    if(ps == 0) ps = 1;
    xb = p * log(y) - lbeta(ps, p) - log(p);
    ans = 0;
    if(xb >= lnsml) {
      ans = exp(xb);
      term = ans * p;
      if(ps != 1) {
        n = fmax(lneps/log(y), 4.0);
        for(i=1 ; i <= n ; i++) {
          xi = i;
          term *= (xi - ps) * y / xi;
          ans += term / (p + xi);
        }
      }
    }
    /* now evaluate the finite sum, maybe. */
    if(q > 1) {
      xb = p * log(y) + q * log1p(-y) - lbeta(p, q) - log(q);
      ib = fmax(xb / lnsml, 0.0);
      term = exp(xb - ib * lnsml);
      c = 1 / (1 - y);
      p1 = q * c / (p + q - 1);
      finsum = 0;
      n = q;
      if(q == n)
        n--;
      for(i=1 ; i<=n ; i++) {
        if(p1 <= 1 && term / eps <= finsum)
          break;
        xi = i;
        term = (q - xi + 1) * c * term / (p + q - xi);
        if(term > 1) {
          ib--;
          term *= sml;
        }
        if(ib == 0) finsum += term;
      }
      ans += finsum;
    }
    if(swap_tail != upper) ans = 1 - ans;
    ans = fmax(fmin(ans, 1.), 0.);
  }
  return ans;
}

double dbinom(long x,long n,double p) {
  double lncomb=lfactorial(n)-lfactorial(x)-lfactorial(n-x);
  return exp(lncomb+x*log(p)+(n-x)*log(1-p));
}

double pbinom(long x,long n,double p) {
  long j;
  double prob=0.0;
  for(j=0;j<=x;j++)
    prob+=dbinom(j,n,p);
  return prob;
}

double dbeta_shift(double x,double *p) {
  double y = p[0];
  double a = p[1];
  double b = p[2];
  return dbeta(x,a,b)-y; 
}

double zeroin(double (*f)(double x,double *params), /* function: f(x,params)        */
              double ax,                            /* lower bound                  */
              double bx,                            /* upper bound                  */
              double *params,                       /* additional parameters for f  */
              double tol,                           /* tolerance on convergence     */
              long maxit) {                         /* maximum number of iterations */
  double a,b,c,fa,fb,fc;
  a = ax;
  b = bx;
  fa = (*f)(a, params);
  fb = (*f)(b, params);
  c = a;
  fc = fa;
  maxit++;
  while(maxit--) {
    double prev_step = b-a;
    double tol_act;
    double p;
    double q;
    double new_step;
    if(fabs(fc) < fabs(fb)) {
      a = b;
      b = c;
      c = a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol_act = 2*DBL_EPSILON*fabs(b) + tol/2;
    new_step = (c-b)/2;
    if(fabs(new_step) <= tol_act || fb == (double)0)
      return b;
    if(fabs(prev_step) >= tol_act && fabs(fa) > fabs(fb) ) {
      register double t1,cb,t2;
      cb = c-b;
      if( a==c ) {
        t1 = fb/fa;
        p = cb*t1;
        q = 1.0 - t1;
      } else {
        q = fa/fc;
        t1 = fb/fc;
        t2 = fb/fa;
        p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
        q = (q-1.0) * (t1-1.0) * (t2-1.0);
      }
      if(p > 0) q = -q;
      else      p = -p;
      if(p < (0.75*cb*q-fabs(tol_act*q)/2) && p < fabs(prev_step*q/2))
        new_step = p/q;
    }
    if(fabs(new_step) < tol_act)
      new_step = (new_step > 0) ? tol_act:-tol_act;
    a = b;
    fa = fb;
    b += new_step;
    fb = (*f)(b, params);
    if((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
      c = a;
      fc = fa;
    }
  }
  /* failed! */
  return b;
}

void binom_bayes(long* x,long* n,
                 double* a,double* b,
                 double* alpha,
                 double* lcl,double* ucl,
                 long* len,
                 long* maxit,double* tol) {
  long j,first,down;
  double lcl_x,ucl_x,lcl_y,ucl_y;
  double y1,y2,y3;
  double px1,px2,sig;
  double mode,xx;
  double x1,x2;
  double lx1,lx2,ux1,ux2;
  double p[3];
  for(j=0;j<len[0];j++) {
    lcl_x=lcl[j];
    ucl_x=ucl[j];
    lcl_y=dbeta(lcl_x,a[j],b[j]);
    ucl_y=dbeta(ucl_x,a[j],b[j]);
    y3=fmax(lcl_y,ucl_y);
    y1=0;
    mode=(a[j]-1)/(a[j]+b[j]-2);
    first=(lcl_y>ucl_y?0:1);
    x1=first?mode:0;
    x2=first?1:mode;
    p[0]=y3; p[1]=a[j]; p[2]=b[j];
    xx=zeroin(dbeta_shift,x1,x2,p,tol[0],maxit[0]);
    if(first) ucl_x= xx;
    else      lcl_x =xx;
    px1=pbeta(lcl_x,a[j],b[j],0);
    px2=pbeta(ucl_x,a[j],b[j],0);
    sig = 1 - px2 + px1;
    down=0;
    while(fabs(sig - 2 * alpha[j]) > tol[0]) {
      y2=(y1+y3)*.5;
      if(down) {
        if(dbeta(lcl_x, a[j], b[j]) < y2) lcl_x = mode;
        lx1=0; lx2=lcl_x;
        if(dbeta(ucl_x, a[j], b[j]) < y2) ucl_x = mode;
        ux1=ucl_x; ux2=1;
      } else {
        if(dbeta(lcl_x, a[j], b[j]) > y2) lcl_x = 0;
        lx1=lcl_x; lx2=mode;
        if(dbeta(ucl_x, a[j], b[j]) > y2) ucl_x = 1;
        ux1=mode; ux2=ucl_x;
      }
      p[0]=y2;
      lcl_x=zeroin(dbeta_shift,lx1,lx2,p,tol[0],maxit[0]);
      ucl_x=zeroin(dbeta_shift,ux1,ux2,p,tol[0],maxit[0]);
      px1=pbeta(lcl_x,a[j],b[j],1);
      px2=pbeta(ucl_x,a[j],b[j],1);
      sig = 1 - px2 + px1;
      if(sig > 2 * alpha[j]) {
        down = 0;
        y3 = y2;
      } else {
        down = 1;
        y1 = y2;
      }
    }
    lcl[j] = lcl_x;
    ucl[j] = ucl_x;
  }
}
