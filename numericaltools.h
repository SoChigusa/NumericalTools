
/**
 * @file numericaltools.h
 * @brief Tools frequently used for numerical analysis
 * @author So Chigusa
 * @date 2018/02/09
 */

#ifndef NUMERICALTOOLS_H
#define NUMERICALTOOLS_H

/* #define _NTOOLS_USE_GETLINE */

#define _NTOOLS_USE_INTEGRATION // and interpolation
/* #define _NTOOLS_CHECK_ASCEND */
/* #define _NTOOLS_USE_PROGRESSBAR */

/*
 * NEED -lMinuit when make
 */
/* #define _NTOOLS_USE_TMINUIT */

/*
 * NEED -lMinuit2 when make
 */
/* #define _NTOOLS_USE_MINUIT2 */

#include <math.h>
#include <iostream>
#include <vector>

#ifdef _NTOOLS_USE_GETLINE
#include <sstream>
#include <string>
#endif

#ifdef _NTOOLS_CHECK_ASCEND
#include <algorithm>
#endif

#ifdef _NTOOLS_USE_TMINUIT
#include <TMinuit.h>
#endif

#ifdef _NTOOLS_USE_MINUIT2
#include <Math/Functor.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <TRandom3.h>
#endif

/**
 * @brief Summarizing class of all functions
 */
class NTools {
public:

  //--------- getline usable for string delimiter ---------
  
#ifdef _NTOOLS_USE_GETLINE
  static void split(const std::string &buffer, const std::string &delim_old,
		    std::vector<std::string> &v, const char delim = ' ');
  static void split(const std::string &buffer, const char delim,
		    std::vector<std::string> &v);
  static void split(const std::string &buffer, const char delim,
		    std::vector<double> &v);
#endif

#ifdef _NTOOLS_USE_INTEGRATION

  //--------- Simpson Integrator (2nd order) ---------
  /**
   * @brief Class for numerical integration using Simpson
   * @details 
   * This is derived by the 2nd order polynomial approximation
   * using successive three discrete points.
   * Please define _NTOOLS_USE_INTEGRATION to use this class.
   */
  class SimpsonIntegrator {
  public:
    static double integrate(std::function<double(double)>, double, double, double);
    static double integrate(std::string, std::vector<double>, std::vector<double>);
  };

  //--------- Boole Integrator (4th order) ---------
  /**
   * @brief Class for numerical integration using Boole's rule
   * @details 
   * This is derived by the 4th order polynomial approximation
   * using successive five discrete points.
   * Please define _NTOOLS_USE_INTEGRATION to use this class.
   */
  class BooleIntegrator {
  public:
    static double integrate(double, std::vector<double>);
  };
  
  //--------- Spline Interpolator ---------
  /**
   * @brief Class for spline interpolation of discrete data
   * @details
   * Please define _NTOOLS_USE_INTEGRATION to use this class.
   */
  class SplineInterpolator {
  private:
    std::vector<double> endpoint;
    std::vector<double> p, q, r, s;
  public:
    void interpolate(std::vector<double>, std::vector<double>);
    double at(double);
    double integrate();
    double integrate(double, double);
    /**
     * @brief Check if interpolation has been done
     * @details 
     * It is recommended to check if this method returns true
     * before using at(double), integrate(), or integrate(double,double) method.
     * @returns true if interpolation has been performed
     */
    bool IsUsed(void) { return !endpoint.empty(); };
  };
#endif
  
  //--------- Progress Bar ---------

#ifdef _NTOOLS_USE_PROGRESSBAR
  class ProgressBar {
  private:
    int pr;
    int end;
  public:
    ProgressBar(int);
    void tick();
  };
#endif
  
  //--------- Wrapper for ROOT Minuit ---------

#ifdef _NTOOLS_USE_TMINUIT
  class Minimization {
  public:
    static void minimize(double &, std::vector<double> & arg_par, std::vector<double> & arg_parerr,
			 void(*arg_fcn)(Int_t &, Double_t *, Double_t &f, Double_t *, Int_t),
			 const std::vector<double> arg_vstart, const std::vector<double> arg_step);
  };

  /* !!!!!!!!! example of the input function !!!!!!!!! */
  /* void chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){ */
  /*   f = par[0]*par[0] + par[0] + 1; */
  /* } */
#endif

#ifdef _NTOOLS_USE_MINUIT2
  class Minimization {
  public:
    static void minimize(double &, std::vector<double> &, std::vector<double> &,
			 const std::function<double(const double *)> &, const std::vector<double> &,
			 const std::vector<double> &, const double,
			 const ROOT::Minuit2::EMinimizerType arg_type = ROOT::Minuit2::kMigrad);
    static void minimize_initial_scan
      (double &, std::vector<double> &, std::vector<double> &,
       const std::function<double(const double *)> &, const std::vector<double> &,
       const std::vector<double> &, const double, const unsigned int,
       const ROOT::Minuit2::EMinimizerType arg_type = ROOT::Minuit2::kMigrad);
  };
#endif
  
};

//--------- Implementation ---------

#ifdef _NTOOLS_USE_GETLINE
void NTools::split(const std::string &buffer, const std::string &delim_old,
		   std::vector<std::string> &v, const char delim) {
  // replacement from delim_old -> delim
  const std::string str_delim{delim};
  std::string s(buffer);
  std::string::size_type pos = 0;
  while ((pos = s.find(delim_old, pos)) != std::string::npos) {
    s.replace(pos, delim_old.length(), str_delim);
    pos += delim_old.length();
  }

  split(s, delim, v);
}

void NTools::split(const std::string &buffer, const char delim,
	     std::vector<std::string> &v) {
  v.resize(0);
  std::stringstream ss(buffer);
  std::string mbuf;
  while(getline(ss, mbuf, delim)) { if(mbuf != "") v.push_back(mbuf); }
}

void NTools::split(const std::string &buffer, const char delim,
	     std::vector<double> &v) {
  v.resize(0);
  std::stringstream ss(buffer);
  std::string mbuf;
  while(getline(ss, mbuf, delim)) { if(mbuf != "") v.push_back(stod(mbuf)); }
}
#endif

#ifdef _NTOOLS_USE_INTEGRATION
/**
 * @brief Numerically integrate an input analytic function
 * @param[in] arg_f input function
 * @param[in] arg_x0 minimum integration range
 * @param[in] arg_x1 maximum integration range
 * @param[in] arg_dx integration step
 * @return double integration result
 */
double NTools::SimpsonIntegrator::integrate(std::function<double(double)> arg_f,
					    double arg_x0, double arg_x1, double arg_dx)
{
  double res = 0.;
  double xmd;
  double xlg = arg_x0;
  double fsm, fmd;
  double flg = arg_f(arg_x0);
  while(xlg < arg_x1) {
    xmd = xlg + 0.5 * arg_dx;
    xlg = xlg + arg_dx;
    fsm = flg;
    fmd = arg_f(xmd);
    flg = arg_f(xlg);
    res += (fsm+4.*fmd+flg)*arg_dx/6.;
  }
  return res;
}

/**
 * @brief Numerically integrate an input discrete data
 * @param[in] arg_flag
 * flag to control integration method:
 * "Simpson" use the Simpson integration method.
 * This method can only be used when equally spaced data is provided.
 * "Spline" first interpolates the date using spline
 * and perform the integration using interpolation parameters.
 * @param[in] arg_x discrete coordinate set
 * @param[in] arg_f discrete data set
 * @return double integration result
 */
double NTools::SimpsonIntegrator::integrate(std::string arg_flag, std::vector<double> arg_x, std::vector<double> arg_f)
{
  const int nmesh = arg_x.size();
  if(arg_f.size() != nmesh) {
    throw "Mismatch between # of mesh of x and y\n";
  }

  double res = 0.;
  if(arg_flag == "Simpson") {
    if(nmesh % 2 == 0) {
      throw "Odd # of date including start and end point is needed\n";
    } else if(nmesh < 3) {
      throw "Too small # of mesh\n";
    }

    const double dx = arg_x[2] - arg_x[0];
    double xmd = arg_x[1];
    for(int i = 1; i < nmesh-1; i=i+2) {
      if(fabs(arg_x[i-1]-(xmd-0.5*dx)) > 1.e-10 ||
	 fabs(arg_x[i]-xmd) > 1.e-10) {
	throw "Non-universal choice of dx\n";
      }
      res += (arg_f[i-1]+4.*arg_f[i]+arg_f[i+1])*dx/6.;
      xmd = xmd+dx;
    }
  }

  else if(arg_flag == "Spline") {
    NTools::SplineInterpolator mySpline;
    mySpline.interpolate(arg_x, arg_f);
    res = mySpline.integrate();
  }

  else {
    throw "Unsupported flag for SimpsonIntegrator\n";
  }
  return res;
}

/**
 * @brief Numerically integrate an input discrete data using Boole's rule
 * @details This method uses 4th order polynomial approximation for the integration.
 * Also, this can only be used when equally spaced data is provided.
 * @param[in] arg_dx data spacing
 * @param[in] arg_f discrete data set
 * @return double integration result
 */
double NTools::BooleIntegrator::integrate(double arg_dx, std::vector<double> arg_f) {
  double res = 0.;
  const int nmesh = arg_f.size();
  static const double f29o24 = 29./24.;
  static const double f71o24 = 71./24.;
  static const double f35o24 = 35./24.;
  static const double to45 = 2./45.;
  res += f29o24*arg_f[0]+f71o24*arg_f[1]+f35o24*arg_f[2]+0.375*arg_f[3]; // 1st, 2nd, 3rd order
  res += 0.375*arg_f[nmesh-4]+f35o24*arg_f[nmesh-3]+f71o24*arg_f[nmesh-2]+f29o24*arg_f[nmesh-1];
  for(int i = 2; i < nmesh-2; ++i) {// 1st order
    res += to45*(7.*arg_f[i-2] + 32.*arg_f[i-1] + 12.*arg_f[i] + 32.*arg_f[i+1] + 7.*arg_f[i+2]);
  }
  res *= 0.25 * arg_dx; // treat over-counting & spacing
  return res;
}

/**
 * @brief Perform spline interpolation of discrete data
 * @details Results of the interpolation are stored in private members
 * and can be accessed through the NTools::SplineInterpolator::at method.
 * @param[in] arg_x discrete coordinate set
 * @param[in] arg_y discrete data set
 */
void NTools::SplineInterpolator::interpolate(std::vector<double> arg_x, std::vector<double> arg_y)
{
  const int n = arg_x.size()-1;
  if(arg_y.size() != n+1) {
    throw "Mismatch between # of x and y points\n";
  } 
  
#ifdef _NTOOLS_CHECK_ASCEND
  else {
    std::vector<double> check;
    std::copy(arg_x.begin(), arg_x.end(), std::back_inserter(check));
    std::sort(check.begin(), check.end());
    if(!std::equal(arg_x.begin(), arg_x.end(), check.begin())) {
      throw "Data should be ascending order in x\n";
    }
  }
#endif

  endpoint.resize(arg_x.size());
  std::copy(arg_x.begin(), arg_x.end(), endpoint.begin());
  std::vector<double> h(n), b(n), d(n), g(n), u(n);
  p.resize(n);
  q.resize(n);
  r.resize(n+1);
  s.resize(n);

  int i;
  for(i=0; i<n; i++) { h[i]=arg_x[i+1]-arg_x[i]; }
  for(i=1; i<n; i++) {
    b[i]=2.*(h[i]+h[i-1]);
    d[i]=3.*((arg_y[i+1]-arg_y[i])/h[i] - (arg_y[i]-arg_y[i-1])/h[i-1]);
  }
  g[1]=h[1]/b[1];
  for(i=2; i<n-1; i++) { g[i]=h[i]/(b[i]-h[i-1]*g[i-1]); }
  u[1]=d[1]/b[1];
  for(i=2; i<n; i++) { u[i]=(d[i]-h[i-1]*u[i-1])/(b[i]-h[i-1]*g[i-1]); }
  r[0]=0; r[n]=0;
  r[n-1]=u[n-1];
  for(i=n-2; i>=1; i--) { r[i]=u[i]-g[i]*r[i+1]; }
  for(i=0; i<n; i++) {
    p[i]=arg_y[i];
    q[i]=(arg_y[i+1]-arg_y[i])/h[i]-h[i]*(r[i+1]+2.*r[i])/3.;
    s[i]=(r[i+1]-r[i])/(3.*h[i]);
  }
}

/**
 * @brief Access the interpolation result
 * @param[in] arg_x coordinate to be accessed
 * @returns data at the point
 */
double NTools::SplineInterpolator::at(double arg_x)
{
  const int n = endpoint.size()-1;
  if(arg_x < endpoint.front() || arg_x > endpoint.back()) {
    std::cout << endpoint.front() << " " << arg_x << " " << endpoint.back() << std::endl;
    throw "Point is not between the interpolation region in SplineInterpolator::at\n";
  }

  int nregion;
  if(arg_x == endpoint[n]) { nregion = n-1; }
  else {
    for(int i = 0; i < n; i++) {
      if(endpoint[i] <= arg_x && arg_x < endpoint[i+1]) { nregion = i; }
    }
  }
  const double xx = arg_x - endpoint[nregion];
  return p[nregion]+q[nregion]*xx+r[nregion]*xx*xx+s[nregion]*xx*xx*xx;
}

/**
 * @brief Integrate the whole range of interpolation
 * @returns integration result
 */
double NTools::SplineInterpolator::integrate()
{
  return integrate(endpoint.front(), endpoint.back());
}

/**
 * @brief Integrate over the given range within the interpolation range
 * @param[in] arg_x0 minimum integration range
 * @param[in] arg_x1 maximum integration range
 * @returns integration result
 */
double NTools::SplineInterpolator::integrate(double arg_x0, double arg_x1)
{
  if(arg_x0 < endpoint.front() || arg_x0 > endpoint.back() ||
     arg_x1 < endpoint.front() || arg_x1 > endpoint.back()) {
    throw "Point is not between the interpolation region in SplineInterpolator::integrate\n";
  } else if(arg_x0 > arg_x1) {
    throw "Integration with x0 > x1 is not supported\n";
  }

  bool fsum = 0;
  double res = 0.;
  double sub_x0, sub_x1;
  const int n = endpoint.size()-1;
  for(int i = 0; i < n; i++) {
    if(fsum == 0) {
      if(endpoint[i+1] > arg_x0) { fsum = 1; }
    }
    if(fsum == 1) {
      sub_x0 = std::max(endpoint[i], arg_x0);
      sub_x1 = std::min(endpoint[i+1], arg_x1);
      auto F = [=](double arg_x){
        double xx = arg_x - endpoint[i];
	return p[i]*xx+q[i]*xx*xx/2.+r[i]*xx*xx*xx/3.+s[i]*xx*xx*xx*xx/4.;
      };
      res += F(sub_x1)-F(sub_x0);
      if(endpoint[i+1] >= arg_x1) { fsum = 2; }
    }
  }
  return res;
}
#endif

#ifdef _NTOOLS_USE_PROGRESSBAR
NTools::ProgressBar::ProgressBar(int arg_e)
: pr(0), end(arg_e)
{
  std::cout << "------------------ ProgressBar -------------------" << std::endl;
}

void NTools::ProgressBar::tick()
{
  pr++;
  if(floor(50.*(pr-1)/end) < floor(50.*pr/end)) { std::cout << "*" << std::flush; }
  if(pr == end) {
    std::cout << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
  }
}
#endif
		 
#ifdef _NTOOLS_USE_TMINUIT
void NTools::Minimization::minimize(double &res_min,
				    std::vector<double> & arg_par, std::vector<double> & arg_parerr,
				    void(*arg_fcn)(Int_t &, Double_t *, Double_t &f, Double_t *, Int_t),
				    const std::vector<double> arg_vstart, const std::vector<double> arg_step)
{
  const int nparam = arg_vstart.size();
  if(nparam != arg_step.size()) {
    throw "Mismatch between start values and steps sizes\n";
  }
  
  TMinuit *gMinuit = new TMinuit(nparam);
  gMinuit->SetPrintLevel(-1);
  gMinuit->SetFCN(arg_fcn);
  
  Int_t ierflg = 0;
  for(int i = 0; i < nparam; i++)
    /* gMinuit->mnparm(i, "name", arg_vstart[i], arg_step[i], 0, 0, ierflg); */
    gMinuit->mnparm(i, "name", arg_vstart[i], arg_step[i], -arg_step[i]*10, arg_step[i]*10, ierflg);

  double par, parerr;
  for(int i = nparam-6; i < nparam-1; i++) {
    gMinuit->GetParameter(i, par, parerr);
    std::cout << par << " ";
  }
  std::cout << std::endl;

  // adjust the error definition for likelihood fits
  double arg_list[2];
  arg_list[0] = 0.5;
  gMinuit->mnexcm("SET ERR", arg_list, 1, ierflg);
  
  arg_list[0] = 500.; // Max Iteration
  arg_list[1] = 1.;   // Tolerance
  gMinuit->mnexcm("SIMplex", arg_list, 2, ierflg);
  for(int i = nparam-6; i < nparam-1; i++) {
    gMinuit->GetParameter(i, par, parerr);
    std::cout << par << " ";
  }
  std::cout << std::endl;
  
  /* gMinuit->Migrad(); */
  arg_par.resize(nparam);
  arg_parerr.resize(nparam);
  for(int i = 0; i < nparam; i++) {
    gMinuit->GetParameter(i, par, parerr);
    arg_par[i] = par;
    arg_parerr[i] = parerr;
  }

  Int_t tmp, tmp2;
  Double_t tmpd;
  arg_fcn(tmp, &tmpd, res_min, arg_par.data(), tmp2);
  
  delete gMinuit;
}
#endif

#ifdef _NTOOLS_USE_MINUIT2
void NTools::Minimization::minimize(double & res_min,
				    std::vector<double> & arg_par,
				    std::vector<double> & arg_parerr,
				    const std::function<double(const double *)> & arg_func,
				    const std::vector<double> & arg_vstart,
				    const std::vector<double> & arg_step,
				    const double arg_tol,
				    const ROOT::Minuit2::EMinimizerType arg_type) {
  size_t nArgs = arg_vstart.size();
  ROOT::Math::Functor F([ &arg_func, &nArgs ](const double* p_phiI) {
      return arg_func(p_phiI);
    }, nArgs);
  ROOT::Minuit2::Minuit2Minimizer min(arg_type);
  /* min.SetMaxFunctionCalls(1000000); */
  min.SetTolerance(arg_tol);
  min.SetFunction(F);

  double step;
  for(int i = 0; i < nArgs; i++) {
    step = arg_step.at(i);
    if(step == 0.) min.SetFixedVariable(i, std::to_string(i), arg_vstart.at(i));
    else min.SetVariable(i, std::to_string(i), arg_vstart.at(i), arg_step.at(i));
  }
  min.Minimize();

  const double *min_x = min.X();
  const double *min_grad = min.Errors();
  arg_par.resize(nArgs);
  arg_parerr.resize(nArgs);
  for(int i = 0; i < nArgs; i++) {
    arg_par[i] = *(min_x + i);
    arg_parerr[i] = *(min_grad + i);
  }
  res_min = min.MinValue();

  // min.PrintResults();

  /* unsigned int nstep = 100; */
  /* double x[100], y[100]; */
  /* min.Scan(1, nstep, x, y, 0., 0.); */
  /* min.Scan(2, nstep, x, y, 0., 0.); */
  /* min.Scan(3, nstep, x, y, 0., 0.); */
  /* min.Scan(4, nstep, x, y, 0., 0.); */
  /* min.Scan(5, nstep, x, y, 0., 0.); */
}
 
void NTools::Minimization::minimize_initial_scan
(double & res_min,
 std::vector<double> & arg_par,
 std::vector<double> & arg_parerr,
 const std::function<double(const double *)> & arg_func,
 const std::vector<double> & arg_vstart,
 const std::vector<double> & arg_step,
 const double arg_tol,
 const unsigned int nscan,
 const ROOT::Minuit2::EMinimizerType arg_type) {
  res_min = 1.e5;
  TRandom3 rgen(1);
  double res_tmp;
  for(int i = 0; i < nscan; i++) {
    std::vector<double> par, parerr;
    std::vector<double> vstart(arg_vstart.size());
    for(int j = 0; j < vstart.size(); j++) {
      vstart[j] = rgen.Gaus(arg_vstart[j], 2.*arg_step[j]);
      // vstart[j] = arg_vstart[j] + 10. * (1. - 2.*rgen.Rndm()) * arg_step[j];
    }
    NTools::Minimization::minimize(res_tmp, par, parerr, arg_func,
				   vstart, arg_step, arg_tol, arg_type);
    if(res_min > res_tmp) {
      res_min = res_tmp;
      arg_par = par;
      arg_parerr = parerr;
    }
  }
  
  /* std::cout << "min = " << res_min << std::endl; */
  /* for(int i = 0; i < arg_par.size(); ++i) { */
  /*   if(arg_step[i] != 0.) { */
  /*     std::cout << i << " = " */
  /* 		<< arg_par[i] << " +- " */
  /* 		<< arg_parerr[i] << std::endl; */
  /*   } */
  /* } */
}

// instance of template
// template void NTools::Minimization::minimize<double (*)(const double *)>
// (double & res_min,
//  std::vector<double> & arg_par,
//  std::vector<double> & arg_parerr,
//  double (*arg_func)(const double *),
//  const std::vector<double> & arg_vstart,
//  const std::vector<double> & arg_step,
//  const double arg_tol,
//  const ROOT::Minuit2::EMinimizerType arg_type);

// template void NTools::Minimization::minimize_initial_scan<double (*)(const double *)>
// (double & res_min,
//  std::vector<double> & arg_par,
//  std::vector<double> & arg_parerr,
//  double (*arg_func)(const double *),
//  const std::vector<double> & arg_vstart,
//  const std::vector<double> & arg_step,
//  const double arg_tol,
//  const unsigned int nscan,
//  const ROOT::Minuit2::EMinimizerType arg_type);
#endif

#endif // NUMERICALTOOLS_H
