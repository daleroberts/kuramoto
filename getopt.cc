#include <unistd.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

template<typename T>
vector<T> split(const T &str, const T &delimiters) {
  vector<T> v;
  typename T::size_type start = 0;
  auto pos = str.find_first_of(delimiters, start);
  while(pos != T::npos) {
    if(pos != start)
      v.emplace_back(str, start, pos - start);
    start = pos + 1;
    pos = str.find_first_of(delimiters, start);
  }
  if(start < str.length())
    v.emplace_back(str, start, str.length() - start);
  return v;
}

int main(int argc, char *argv[]) {
  string graphfile;
  string parameters;
  extern char *optarg;
  extern int optind, optopt;
  int c;
  while ((c = getopt(argc, argv, "p:g:")) != -1) {
    switch(c) {
      case 'p':
        char *key, *value, *str;
        str = strdup(optarg);
        while (true) {
          key = strsep(&str, "=");
          value = strsep(&str, ",");
          if ((key == NULL) || (value == NULL))
            break;
          printf("key %s value %s\n", key, value);
        }
        free(str);
        break;
      case 'g':
        graphfile.assign(optarg);
        break;
      case '?':
        cerr << "Unrecognized option: -"  << optopt << endl;
      default:
        cout << endl;
        abort();
    }
  }

  cout << "graph file: " << graphfile << endl;

  return 0;
}
