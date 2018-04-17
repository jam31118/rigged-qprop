#ifndef PARAMETER_H
#define PARAMETER_H

#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::istringstream;
using std::stringstream;

///Simple Object representing a  parameter.
struct parameter {
  parameter(string n, string t, string v)
    : _name(n), _type(t), _value(v)
  {};
  string _name, _type, _value;
};

///STL type functor for compariosn of parmeter name
struct compareNameWith : public std::unary_function< parameter, bool > {
  compareNameWith(string n) 
    : _name(n)
  {};
  string _name;
  bool  operator() (parameter cmp) { return (cmp._name == _name); };
};

///Class for reading parameters from a file.
class parameterListe {
  std::vector<parameter> _param_list;
public:
  parameterListe();   
  parameterListe(string file_name="parameters");
  void setDouble(string name, double value);
  double getDouble(string name) const;
  string getString(string name) const;
  long getLong(string name) const;
  bool getBool(string name) const;
};

#endif
