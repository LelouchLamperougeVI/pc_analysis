/* csvParse.cpp
 * Takes a csv file, parse through it with the specified format identifier and save the results as a .mat file for use in MATLAB
 *
 * Compilation: g++ -lpthread csvParse.cpp -o csvParse
 * Compiled with g++ 9.1.0 and tested on Arch Linux kernel 5.2.8
 *
 *
 * By HaoRan Chang, Ph.D. candidate
 * Canadian Centre for Behavioural Neuroscience,
 * Department of Neuroscience, University of Lethbridge,
 * Lethbridge, Alberta, Canada
 *
 */

#include <iostream>
#include <cctype>
#include <cstring>
#include <fstream>
#include <string>
#include <thread>
#include <mutex>

#define CHAR_LEN 80

namespace CSV {
  std::ifstream stream;
  char format[CHAR_LEN];
  std::mutex mute;
  size_t d_size; //the size of each datum on memory
  int d_size; //the number of elements in each line

  class parser;
  class data;
};

class CSV::data {
  public:
    void *datum;

    data() {}

    data(int size) {
      datum = malloc(size);
    }
};

class CSV::parser {
  public:
    char buff[CHAR_LEN];

    void parsel() { //parse line

    }

    void print() {
      printf("%s\n", buff);
    }

    parser() {
      int num_d = 0, num_s = 0;
      char *id = format;

      do {
        switch(*id) {
          case 'd':
            num_d++;
          case 's':
            num_s++;
        }
      } while(*id++ != '\0');

      d_size = num_d*sizeof(double) + num_s*CHAR_LEN*sizeof(char);
      d_len = num_d + num_s;
    }

  //private:
    void readline() {
      mute.lock();
      stream.getline(buff, CHAR_LEN);
      mute.unlock();
    }
};

int main () {
  CSV::stream.open("test.csv", std::ifstream::in);
  CSV::parser temp;
  temp.readline();
  temp.print();
  CSV::parser temp2;
  temp2.readline();
  temp2.print();
  CSV::stream.close();

  CSV::data test(sizeof(double) + CHAR_LEN*sizeof(char));
  void *pt = test.datum;
  *(double *) test.datum = 69;
  test.datum = static_cast<double *>(test.datum) + 1;
  strcpy( (char *) test.datum, "lolz");
  std::cout << *(double *) pt << std::endl;
  pt = static_cast<double *>(pt) + 1;
  std::cout << (char *) pt << std::endl;
}
