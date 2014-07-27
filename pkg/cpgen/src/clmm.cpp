/*
// clmm.cpp
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

#include "clmm.h"
#include "printer.h"

typedef vector<MCMC<base_methods_st> > mcmc_st;
typedef vector<MCMC<base_methods_mp> > mcmc_mp;


SEXP clmm(SEXP yR, SEXP XR, SEXP par_XR, SEXP list_of_design_matricesR, SEXP par_design_matricesR, SEXP par_mcmcR, SEXP verboseR, SEXP threadsR){

int threads = as<int>(threadsR); 
int verbose = as<int>(verboseR); 


//omp_set_dynamic(0);
omp_set_num_threads(threads);

Eigen::setNbThreads(1);
Eigen::initParallel();

Rcpp::List list_of_phenotypes(yR);
int p = list_of_phenotypes.size();
int index;

printer prog(p / threads);

mcmc_st vec_mcmc_st;
mcmc_mp vec_mcmc_mp;


if((p>1) | (threads==1)) {

  vec_mcmc_st.resize(p);
  for(mcmc_st::iterator it = vec_mcmc_st.begin(); it != vec_mcmc_st.end(); it++) {

    index = it - vec_mcmc_st.begin();
    it->populate(list_of_phenotypes[index], XR, par_XR, list_of_design_matricesR ,par_design_matricesR, par_mcmcR,index);
    it->initialize();

  }

  if ((p > 1) & verbose) { prog.initialize(); }


// this looks easy - the work was to allow this step to be parallelized
#pragma omp parallel for 
  for(unsigned int i=0;i<vec_mcmc_st.size();i++){

    vec_mcmc_st.at(i).gibbs();

// verbose

    if(p>1) {

      if(verbose) { 
   
        if(omp_get_thread_num()==0) {

          prog.DoProgress(); 

        }

      }

    }

  }


  Rcpp::List Summary;
  for(mcmc_st::iterator it = vec_mcmc_st.begin(); it != vec_mcmc_st.end(); it++) {

    it->summary();
    Summary[it->get_name()] = it->get_summary();     

  }

  return Summary;

} else {

// if the number of threads is larger than and the number of phenotypes is equal to 1 
// the function runs one parallelized Gibbs Sampler - Fernando et al. 2014
    vec_mcmc_mp.resize(p);
    for(mcmc_mp::iterator it = vec_mcmc_mp.begin(); it != vec_mcmc_mp.end(); it++) {

      index = it - vec_mcmc_mp.begin();
      it->populate(list_of_phenotypes[index], XR, par_XR, list_of_design_matricesR ,par_design_matricesR, par_mcmcR,index);
      it->initialize();
      it->gibbs();

    }


    Rcpp::List Summary;
    for(mcmc_mp::iterator it = vec_mcmc_mp.begin(); it != vec_mcmc_mp.end(); it++) {

      it->summary();
      Summary[it->get_name()] = it->get_summary();     

    }


    return Summary;

  }


}



