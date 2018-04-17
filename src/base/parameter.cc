#include "parameter.hh"

///reads from file "parameters"
parameterListe::parameterListe() {
  parameterListe("parameters");
};
// parameterListe::parameterListe() {
//   std::ifstream input("parameters");
//   if(!input)      {
//     std::cerr << "file not read in parameterListe::parameterListe()" << "\n";     
//     exit(-1);
//   };
//   char line[80];
//   while(input.getline(line,80)) {
//     if (line[0] != '#') {
//       string strline(line);    
//       std::istringstream  in_stream(strline);
//       string name, type, value;
//       in_stream >> name >> type >> value;
//       _param_list.push_back(parameter(name, type, value));
//     };
//   };
        
// };

///reads from file file_name
parameterListe::parameterListe(string file_name) {
  std::ifstream input(file_name.c_str());
  if(!input)  {
    std::cerr << "file not read in parameterListe::parameterListe(string file_name)" << "\n";     
    exit(-1);
  };
  char line[1024];
  while(input.getline(line, 1024)) {
    // strip comments from line
    for (size_t i(0); line[i]!='\0'; i++) {
      if (line[i]=='#') {
	line[i]='\0';
	break;
      };        
    };
    // copy the c_string to a nice string
    string strline(line);
    // check if line has exactly three entries
    size_t num_entries(0);
    istringstream  test_stream(strline);
    string throwaway;
    while (test_stream >> throwaway) {
      num_entries++;
    };
    if (num_entries==3) {
      std::istringstream  in_stream(strline);
      string name, type, value;
      in_stream >> name >> type >> value;
      // find out if name was used more than once
      compareNameWith mycmp(name);
      std::vector<parameter>::const_iterator it;
      it = find_if(_param_list.begin(), _param_list.end(), mycmp); 
      if (it != _param_list.end()) {
	cerr << name << " appears more than once in parameter file. You're gonna have a bad time!" << endl;
      }
      else
	_param_list.push_back(parameter(name, type, value));
      cout << name << " " << type << " " << value << endl; 
    }
    else if (num_entries==0) {
      // do nothing for empty line 
    }
    else {
      cerr << strline << " is not a well formed entry and discarded" << endl;
    };
  };
        
};

///set the parameter name of type double to a new value if it exists and return 255 otherwise
void parameterListe::setDouble(string name, double value) {
  compareNameWith mycmp(name);
  std::vector<parameter>::iterator it;
  it = find_if(_param_list.begin(), _param_list.end(), mycmp); 
  if (it == _param_list.end()) {
    std::cerr << "Paramter " << name << " nicht gefunden in parameterListe::getDouble" << std::endl;
    exit(-1);
  }
  else if (it->_type != "double") {
    std::cerr << "Falscher Typ in parameterListe::getDouble" << std::endl;
    exit(-1);
  }  
  else {
    stringstream bla;
    bla << value; 
    it->_value = bla.str();
  };
};

///return the parameter name of type double if it exists and return 255 otherwise
double parameterListe::getDouble(string name)  const {
  compareNameWith mycmp(name);
  std::vector<parameter>::const_iterator it;
  it = find_if(_param_list.begin(), _param_list.end(), mycmp); 
  if (it == _param_list.end()) {
    std::cerr << "Paramter " << name << " nicht gefunden in parameterListe::getDouble" << std::endl;
    exit(-1);
  }
  else if (it->_type != "double") {
    std::cerr << "Falscher Typ in parameterListe::getDouble" << std::endl;
    exit(-1);
  }  
  else {
    stringstream bla(it->_value);
    double blub;
    bla >> blub;
    return blub;
  };
};

long parameterListe::getLong(string name)  const {
  compareNameWith mycmp(name);
  std::vector<parameter>::const_iterator it;
  it = find_if(_param_list.begin(), _param_list.end(), mycmp); 
  if (it == _param_list.end())
    {
      std::cerr << "Paramter " << name << " nicht gefunden in parameterListe::getLong" << std::endl;
      exit(-1);
    }
  else if (it->_type != "long")
    {
      std::cerr << "Falscher Typ in parameterListe::getLong" << std::endl;
      exit(-1);
    }  
  else
    {
      stringstream bla(it->_value);
      int blub;
      bla >> blub;
      return blub;
    };
};

string parameterListe::getString(string name)  const {
  compareNameWith mycmp(name);
  std::vector<parameter>::const_iterator it;
  it = find_if(_param_list.begin(), _param_list.end(), mycmp); 
  if (it == _param_list.end()) {
    std::cerr << "Paramter " << name << "  nicht gefunden in parameterListe::getString" << std::endl;
    exit(-1);
  }
  else if (it->_type != "string") {
    std::cerr << "Falscher Typ in parameterListe::getString" << std::endl;
    exit(-1);
  }  
  else {
    stringstream bla(it->_value);
    string blub;
    bla >> blub;
    return blub;
  };
};

bool parameterListe::getBool(string name)  const {
  compareNameWith mycmp(name);
  std::vector<parameter>::const_iterator it;
  it = find_if(_param_list.begin(), _param_list.end(), mycmp); 
  if (it == _param_list.end()) {
    std::cerr << "Paramter " << name << "  nicht gefunden in parameterListe::getBool" << std::endl;
    exit(-1);
  }
  else if (it->_type != "bool") {
    std::cerr << "Falscher Typ in parameterListe::getBool" << std::endl;
    exit(-1);
  }  
  else {
    stringstream bla(it->_value);
    bool blub;
    bla >> blub;
    return blub;
  };
};
