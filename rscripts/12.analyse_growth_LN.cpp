#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fExpCpp(NumericVector x, NumericVector param, int t_diff, double max_diam){
  int n = x.size();
  NumericVector q = (max_diam - x) * (1-exp(-t_diff*(param[0]*log(x)+param[1])));
  NumericVector z(n);
  double mean_z = 0;
  int nq = 0;
  
  for(int i=0;i<n;i++){
    if(q[i]>0) {
      z[i] = log(q[i]);
      mean_z = mean_z + log(q[i]);
      nq++;
    } else {
      z[i] = 0;
    }
  }
  
  if(nq != 0){
    
    mean_z = mean_z / nq;
    
    for(int i=0;i<n;i++){
      if(z[i] == 0){
        z[i] = mean_z;
      }
    }
    
  }
  return z;
  
}

// [[Rcpp::export]]
double quadTrapezCpp_1(NumericVector q, double h, int nx){
  double int_q = sum(q) - 0.5 * (q[0]+q[nx-1]);
  return int_q * h;
}

//[[Rcpp::export]]
void my_dlnorm(NumericVector means, NumericVector sds, double h, int nx, NumericMatrix mat){
    int nrow = mat.nrow(), ncol = mat.ncol();
    for( int i=0; i<nrow; i++){
    	for( int j=i; j<ncol; j++){
      		mat(i,j) =  R::dlnorm(h*(j-i),means[i],sds[i],0);
    	}
    }
}


// [[Rcpp::export]]
void InverseFentonWilkinson(NumericVector Mu, NumericVector Var,  NumericVector mx, NumericVector vx) {
    vx = log(10*(exp(Var)-1)+1);
    mx = Mu-(vx-Var)/2-log(10);
}

// [[Rcpp::export]]
NumericMatrix IPMpdfTreeGrowthCpp(NumericVector x, NumericVector y, NumericVector param, int t_diff, double max_diam, int nx, double h){
  NumericMatrix mat(x.size(),x.size());
  NumericVector f_gr_10y = fExpCpp(x,param,t_diff,max_diam);
  NumericVector variance_10y = exp(param[2] + param[3]*log(x));
  NumericVector f_gr(f_gr_10y.size());
  NumericVector variance(variance_10y.size());
  InverseFentonWilkinson(f_gr_10y,variance_10y, f_gr, variance);
  my_dlnorm(f_gr,sqrt(variance), h, nx, mat);
  return mat;
}

// [[Rcpp::export]]
NumericMatrix IPMpdfTreeGrowth10yCpp(NumericVector x, NumericVector y, NumericVector param, int t_diff, double max_diam, int nx, double h){
  NumericMatrix mat(x.size(),x.size());
  NumericVector f_gr = fExpCpp(x,param,t_diff,max_diam);
  NumericVector variance = exp(param[2] + param[3]*log(x));
  my_dlnorm(f_gr,sqrt(variance), h, nx, mat);
  return mat;
}

// [[Rcpp::export]]
NumericVector CoefTimesVarCpp(Rcpp::List param, NumericVector z, NumericVector s = 0, NumericVector spba = 0, std::string type = "survival") {
  
  int n = Rcpp::as<NumericVector>(param[0]).size();
  NumericVector res(n);
  
  for(int i = 0;i<n;i++){
    
    if (type == "ingrowth"){
      
      res[i] = as<NumericVector>(param["intercept"])[i] + 
        as<NumericVector>(param["rain"])[i]*z[0] + 
        as<NumericVector>(param["temper"])[i]*z[1] +
        as<NumericVector>(param["basal.area"])[i]*z[2] +
        as<NumericVector>(param["anom.rain"])[i]*z[3] +
        as<NumericVector>(param["anom.temper"])[i]*z[4] +
        as<NumericVector>(param["rain.anom.temper"])[i]*z[0]*z[4];
      
      if(s.size() > 1){
        res[i] += as<NumericVector>(param["saplings"])[i]*s[i];
      }
      
      res[i] = res[i] * (10000/(M_PI*25));
      
    }else if (type == "growth" | type == "survival"){
      
      res[i] = as<NumericVector>(param["intercept"])[i] + 
        as<NumericVector>(param["rain"])[i]*z[0] + 
        as<NumericVector>(param["temper"])[i]*z[1] +
        as<NumericVector>(param["basal.area"])[i]*z[2] +
        as<NumericVector>(param["anom.rain"])[i]*z[3] +
        as<NumericVector>(param["anom.temper"])[i]*z[4] +
        as<NumericVector>(param["rain.anom.temper"])[i]*z[0]*z[4];
      
    }else if (type == "recruitment"){
      
      res[i] = as<NumericVector>(param["intercept"])[i] + 
        as<NumericVector>(param["rain"])[i]*z[0] + 
        as<NumericVector>(param["temper"])[i]*z[1] +
        as<NumericVector>(param["basal.area"])[i]*z[2] +
        as<NumericVector>(param["anom.rain"])[i]*z[3] +
        as<NumericVector>(param["anom.temper"])[i]*z[4] +
        as<NumericVector>(param["rain.anom.temper"])[i]*z[0]*z[4] + 
        as<NumericVector>(param["basal.area.sp"])[i]*spba[i];
      
    }
  }
  return res;
}
