#include <Rcpp.h>

#include <fstream>
#include <sstream>
#include <string>

using namespace Rcpp;
using namespace std;

//[[Rcpp::export]]
StringVector my_split(string s, char delim)
{
  StringVector result;
  stringstream ss(s);
  string item;
  while (getline(ss, item, delim)) {
    result.push_back(item);
  }
  return result;
}

