#ifndef _Strings_h_
#define _Strings_h_
using namespace std;

/*
 * String parsing functions
 */

// trim from start
static inline std::string &ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
  return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
  return ltrim(rtrim(s));
}

// parse flag-based input ("-a AAA -b BBB")
static inline bool parseInputWithFlag(int argc, char **argv, char flag, string *value) {
  bool bopt = false;

  for(int i=1; i < argc; i++) {
    if(argv[i][0] == '-' && argv[i][1] == flag) {
      if(argv[i+1][0] != '-') {
        if(value != NULL) *value = argv[i+1];
        bopt = true;
        i++;
        break;
      }
    }
  }

  return bopt;
}

//parse next whitespace separated value from buffer
static inline bool parseNextValue(string *buffer, string *value) {
  int i = buffer->find_first_not_of(" \t,");
  if(i) buffer->erase(0, i);
  i = buffer->find_first_of(" \t,\n\r");

  if(i == -1) {
    *value = buffer->substr(0, i);  //Value of string::npos for second argument grabs entire string
    buffer->erase(0, buffer->length());
    return false;
  } else {
    *value = buffer->substr(0, i);
    buffer->erase(0, i+1);
  }

  return true;
}

// 
static inline bool starts_with(string *buffer, string match) {
  int n = match.length();
  string begin = buffer->substr(0, n);
  if(begin == match) {
    return true;
  }
  return false;
}

static inline bool ends_with(string *buffer, string match) {
  int n = match.length();
  string ending = buffer->substr(buffer->length()-n, n);
  if(ending == match) {
    return true;
  }
  return false;
}

#define charToDouble(str) strtod(str, NULL)
#define stringToFloat(cppstr) strtof(cppstr.c_str(), NULL)
#define stringToDouble(cppstr) strtod(cppstr.c_str(), NULL)
#define stringToInt(cppstr) atoi(cppstr.c_str())



#endif
