//#include "openfile.h"
#include <experimental/filesystem>
#include <iostream>

namespace fs = std::experimental::filesystem;
using namespace std;

int main(int argc, char *argv[]) {
  fs::path path (argv[1]);

  string temp = fs::canonical(path);

  temp = temp.substr(0, temp.rfind(".")) + ".mat";
  cout << temp << endl;
}
