/*
// mcmc.h
// Claas Heuer, June 2014
//
// Copyright (C)  2014 Claas Heuer
//
// This file is part of cpgen.
//
// cpgen is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// cpgen is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// in file.path(R.home("share"), "licenses").  If not, see
// <http://www.gnu.org/licenses/>.
*/

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<Eigen::MatrixXd> MapMatrixXd;
typedef Eigen::MappedSparseMatrix<double> MapSparseMatrixXd;
typedef Eigen::Map<Eigen::VectorXd> MapVectorXd;
typedef Eigen::Map<Eigen::ArrayXd> MapArrayXd;

// Taken from: http://gallery.rcpp.org/articles/sparse-iterators/
typedef MapSparseMatrixXd::InnerIterator InIterMat;
typedef Eigen::SparseVector<double> SpVec;
typedef SpVec::InnerIterator InIterVec;

typedef std::vector<std::map<std::string, int> > mp_container;

//#ifndef _CRAN
// #include "mt_sampler.h"
//#else
// #include "R_sampler.h"
//#endif

#include "mt_sampler.h"
#include "effects.h"
#include "printer.h"

////////////////
// MCMC class //
////////////////

template<class F>
class MCMC {

private:

  int niter;
  int iter;
  int  burnin;
  int  p;
  bool full_output;
  bool verbose;
  bool initialized;

  double scale_e;
  double df_e;

  sampler mcmc_sampler;

public:

  Rcpp::List list_of_design_matrices;
  Rcpp::List summary_list;
  double var_y;
  double mu_y;
  double var_beta;
  double post_var_e;
  double mu;
 
  int effiter;
  int n_random;

  int n;

  uint32_t seed;
 
  std::string name;

  VectorXd y;
  VectorXd ycorr;
  VectorXd mean_var_e;

  mp_container thread_vec;

// this obejct determines whethter we run in parallel or not
  base_methods_abstract * my_base_functions;

  vector<effects> model_effects;

  double var_e;

  bool has_na;
  vector<int> isna;

  void populate(SEXP y_from_R, SEXP X_from_R, SEXP par_fixed_from_R, SEXP list_of_design_matrices_from_R, SEXP par_design_matrices_from_R, SEXP par_from_R, int phenotype_number);
  inline void initialize();
  inline void sample_random();
  inline void sample_residual();
  inline void finish_iteration();
  inline void gibbs();
  inline void summary();
  std::string get_name();
  Rcpp::List get_summary();
// taken from http://stackoverflow.com/questions/307082/cleaning-up-an-stl-list-vector-of-pointers
//  ~MCMC() { while(!model_effects.empty()) delete model_effects.back(), model_effects.pop_back() ;};
  MCMC() : my_base_functions(new F) {};
  ~MCMC(){ delete my_base_functions; };

};

template<class F>
void MCMC<F>::populate(SEXP y_from_R, SEXP X_from_R, SEXP par_fixed_from_R, SEXP list_of_design_matrices_from_R, SEXP par_design_matrices_from_R, SEXP par_from_R, int phenotype_number) {

  Rcpp::List par(par_from_R);

  niter = Rcpp::as<int>(par["niter"]);
  burnin = Rcpp::as<int>(par["burnin"]);
  full_output = Rcpp::as<bool>(par["full_output"]);
  verbose = Rcpp::as<bool>(par["verbose"]);
  scale_e = Rcpp::as<double>(par["scale_e"]);
  df_e = Rcpp::as<double>(par["df_e"]);

  seed = Rcpp::as<uint32_t>(par["seed"]);

  std::ostringstream oss; 
  oss << phenotype_number + 1; 
  name = "PHENOTYPE_" + oss.str();

  list_of_design_matrices = Rcpp::List(list_of_design_matrices_from_R);
  MapVectorXd y_temp = MapVectorXd(as<MapVectorXd> (y_from_R));
  y = y_temp;

  mean_var_e = VectorXd::Zero(niter);

  effiter=0;
  post_var_e=0;
  iter=0;

  n = y.size();
  n_random = list_of_design_matrices.size() -1;
  if(n_random<1) n_random=1;

// check for NAs in pehnotype vector
// comparison x!=x yields TRUE if x=NA

  has_na=0;
  mu = 0;

  for(int i=0;i<y.size();i++) { if( y(i)!=y(i) ){ isna.push_back(i);} else { mu+=y(i);} }
//Rcout << endl << "na size: " << isna.size() << endl;
  if(isna.size() > 0) { has_na=1;}

  mu = mu / (y.size() - isna.size());  

// compute variance; missing at random for NAs -- FIXME crucial part
// FIXME assign: if( y(i)!=y(i)  or if( y(i)!=0 ??
  for(int i=0;i<y.size();i++) { if( y(i)!=y(i) ) {y(i) = 0;} else {var_y += (y(i) - mu) * (y(i) - mu);} }
  var_y = var_y / (y.size() - isna.size() - 1);
  var_e = var_y - (var_y / (2 * n_random));

  ycorr = y;

/////////////////////
// Multithreading //
////////////////////

  int n_threads;

// container to store start and length for OpenMp threads - Credit: Hao Cheng (Iowa State University)
// get the number of threads
#pragma omp parallel
{
if(omp_get_thread_num()==0) { n_threads = omp_get_num_threads(); }
}


thread_vec.resize(n_threads);


for(int i=0;i<n_threads;i++) {

  thread_vec.at(i)["start"] = i * n / n_threads;
  if(i==(n_threads-1)){thread_vec.at(i)["end"] = n;} else {
  thread_vec.at(i)["end"] = (i+1) * n / n_threads;}
  thread_vec.at(i)["length"] = thread_vec.at(i)["end"] - thread_vec.at(i)["start"];

}



// populate with effects
// include fixed effect

  Rcpp::List par_fixed(par_fixed_from_R);
  model_effects.push_back(effects(X_from_R, ycorr.data(),par_fixed, &var_e, var_y,n_random,niter, burnin, full_output)); 

// random effects

  Rcpp::List par_design_matrices(par_design_matrices_from_R);

  for(int i=0;i<list_of_design_matrices.size();i++){

    SEXP temp_list_sexp = par_design_matrices[i];
    Rcpp::List temp_list(temp_list_sexp);

    model_effects.push_back(effects(list_of_design_matrices[i],ycorr.data(),temp_list, &var_e, var_y,n_random,niter, burnin, full_output)); 

  }

 initialized=0;

 mcmc_sampler.set_seed(seed+phenotype_number);

//  Rcout << endl << "MCMC I am single_thread" << endl;


}


template<class F>
void MCMC<F>::initialize() {


    for(vector<effects>::iterator it = model_effects.begin(); it != model_effects.end(); it++) {

    it->initialize(my_base_functions);

    }

}


template<class F>
void MCMC<F>::gibbs() {

vector<effects>::iterator it;  

int progress = 0;
printer prog(niter);


  for(int gibbs_iter=0; gibbs_iter<niter; gibbs_iter++){


    for(it = model_effects.begin(); it != model_effects.end(); it++) {

    it->sample_effects(mcmc_sampler,my_base_functions,thread_vec);

    }

    for(it = model_effects.begin(); it != model_effects.end(); it++) {

    it->sample_variance(mcmc_sampler, gibbs_iter);

    }



// sample residual variance

   var_e = (ycorr.matrix().squaredNorm() + scale_e * df_e) / mcmc_sampler.rchisq(n + df_e);
   mean_var_e(gibbs_iter)=var_e;

// residual noise to NAs -- Reference: de los Campos (2009) - BLR
    if(has_na) 
     {
       for (unsigned int i=0;i<isna.size();i++) {ycorr(isna[i]) = mcmc_sampler.rnorm(0,sqrt(var_e));}
     }

// posterior means

    if(gibbs_iter>burnin) {

      for(it = model_effects.begin(); it != model_effects.end(); it++) {

      it->update_means();

      }

    post_var_e += var_e;
    effiter++;

    }

    if (verbose) {

      progress++;
      prog.DoProgress(progress);

    }

//    if(verbose) { Rcout << endl << "Iteration: " << gibbs_iter << "   s2e: " << var_e; }

  }

//  if(verbose) { Rcout << endl; }


}


template<class F>
void MCMC<F>::summary() {


  summary_list["Residual_Variance"] = Rcpp::List::create(Rcpp::Named("Posterior_Mean") = post_var_e / effiter,
			      Rcpp::Named("Posterior") = mean_var_e,
			      Rcpp::Named("scale_prior") = scale_e,
			      Rcpp::Named("df_prior") = df_e);

  VectorXd yhat = VectorXd::Zero(n);
  for(vector<effects>::iterator it = model_effects.begin(); it != model_effects.end(); it++) {

    yhat += it->predict(effiter);

  }
  

  summary_list["Predicted"] = yhat;

  int count = 1;

  std::string list_name;

  for(vector<effects>::iterator it = model_effects.begin(); it != model_effects.end(); it++) {
  
    std::ostringstream oss; 
    oss << count;
    list_name = "Effect_" + oss.str();
    summary_list[list_name] = it->get_summary(effiter);
    count++;
  }
  

}


template<class F>
std::string MCMC<F>::get_name() {


  return name;
  

}


template<class F>
Rcpp::List MCMC<F>::get_summary() {

 return summary_list;

}


