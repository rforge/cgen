#ifndef _PRINTER_H
#define _PRINTER_H
#include <stdio.h>
#include <iomanip> 
#include <Rcpp.h>


// taken from: http://www.codeproject.com/Tips/537904/Console-simple-progress

class printer {

private:

int percent;
int previous_percent;
static const int width = 40;
int pos;
int previous_pos;
int total;
bool initialized;
std::ostringstream oss;
std::string position; 

public:

inline void DoProgress( int step)
{



    percent = ( step * 100 ) / total;

    if(percent>100) { percent=100; }

    if(percent > previous_percent) {

      if(!initialized) { init(); }

      fflush(stdout);
      previous_pos = pos;
      pos = ( step * width ) / total;  
      step = pos - previous_pos;

      for ( int i = 0; i < step; i++ )  Rcpp::Rcout << "=";

      for(int i=0;i<(width-pos+2);i++) { Rcpp::Rcout << "\033[C"; }

      Rcpp::Rcout << percent << "%" << "\033[" << pos+2 << "G";

      previous_percent = percent;

      if(pos == width) { Rcpp::Rcout << std::endl << std::endl; }

   }


};

inline void init() {


//cout << "\033[1G" << n*10 << "\033[5G" << flush;

      Rcpp::Rcout << std::endl << "[";

      //fill progress bar with =
      for ( int i = 0; i < width; i++ )  Rcpp::Rcout << " ";

      Rcpp::Rcout << "]" << " " << "0%"<< "\033[2G" <<  std::flush;

      initialized = 1;

};


printer(int t) : total(t), pos(0), previous_pos(0), previous_percent(0), initialized(0) {} 



};



#endif
