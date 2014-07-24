/*
// mt_sampler.h
// Claas Heuer, June 2014
//
// Copyright (C)  2014 Claas Heuer
//
// This file is part of cgenpp.
//
// cgenpp is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// cgenpp is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// in file.path(R.home("share"), "licenses").  If not, see
// <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <random>

class sampler {

private:
  std::mt19937 gen;

public:
  double rnorm(double mean, double sd){

    return std::normal_distribution<double>(mean,sd)(gen);

  };

  double rchisq(int df){

    return std::chi_squared_distribution<double>(df)(gen);

  };

  void set_seed(uint32_t seed){

  gen.seed(static_cast<uint32_t> (seed));

  };

  void set_seed(){

  std::random_device rd;
  gen = std::mt19937(static_cast<uint32_t> (rd())); 

  };

  void check_sampler(){ std::cout << std::endl << " I am mt19937_C++ " << std::endl;};

};

