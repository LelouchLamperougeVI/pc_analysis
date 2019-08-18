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
 */

#include <iostream>
#include <ctime>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <mutex>

#define CHAR_LEN 80
#define BUFF_LEN 500
#define ARGS_OFFSET 3

namespace CSV {
  std::ifstream stream; //file stream
  char *format; //format identifier
  char **tokens; //search tokens
  std::mutex mute; //mutex guard
  size_t d_size; //the size of each datum on memory
  int d_len;

  void initialize(char *fn, char *iformat, char **itokens);
  void exit();

  class parser;
  class data;
};

void CSV::initialize(char *fn, char *iformat, char **itokens) {
  format = iformat;
  tokens = itokens;
  stream.open(fn, std::ifstream::in);
  int num_d = 0, num_s = 0;
  char *id = format;
  do {
    switch(*id) {
      case 'd':
        num_d++;
        break;
      case 's':
        num_s++;
        break;
    }
  } while(*id++ != '\0');
  d_size = num_d*sizeof(double) + num_s*CHAR_LEN*sizeof(char);
  d_len = num_d + num_s;
}
void CSV::exit() {
  stream.close();
}

class CSV::data {
  public:
    void *datum;

    data() {
      datum = malloc(d_size);
      write_head = datum;
      currPos = format;
    }

    void free_mem() {
      free(datum);
    }

    void push(void *new_data) {
      SKIP: switch(*currPos++) {
        case 'd':
          *(double *) write_head = *(double *) new_data;
          write_head = static_cast<double *>(write_head) + 1;
          break;
        case 's':
          strcpy( (char *) write_head, (char *) new_data );
          write_head = static_cast<char *>(write_head) + CHAR_LEN;
          break;
        case 'x':
          goto SKIP; //so I used 'GOTO', big deal
      }
    }

    void *pop() {
      SKIP: switch(*--currPos) {
        case 'd':
          write_head = static_cast<double *>(write_head) - 1;
          break;
        case 's':
          write_head = static_cast<char *>(write_head) - CHAR_LEN;
          break;
        case 'x':
          goto SKIP;
      }
      return write_head;
    }

  private:
    void *write_head;
    char *currPos;
};

class CSV::parser {
  public:
    char buff[BUFF_LEN];
    std::vector<data> data_vect;
    char *sp = format; //stack pointer, old habits from assembly will never change...

    parser() {
      while(*++sp != '\0') {}
    }

    ~parser() {
      for(int i = 0; i < data_vect.size(); i++)
        data_vect[i].free_mem();
      data_vect.clear();
    }

    int readline() {
      mute.lock();
      stream.getline(buff, BUFF_LEN);
      if(stream.eof() || stream.bad() || stream.fail()) {
        mute.unlock();
        return 1;
      }
      mute.unlock();
      parsel(buff);
      return 0;
    }

    void *consolidate() { //pop out column from last data element
        SKIP: switch(*--sp) {
        case 'd': {
          double *column = new double[data_vect.size()];
          for(int i = 0; i < data_vect.size(); i++) {
            double *temp = (double *) data_vect[i].pop();
            column[i] = *temp;
          }
          return column;
        }
        case 's': {
          char **column = new char*[data_vect.size()];
          for(int i = 0; i < data_vect.size(); i++) {
            column[i] = new char[CHAR_LEN];
            column[i] = (char *) data_vect[i].pop();
          }
          return column;
        }
        case 'x':
          goto SKIP;
        default:
          return NULL;
      }
    }

  private:
    void parsel(char *line) { //parse line
      int pos = 0;
      double d;
      data temp;
      char *tok = strtok(line, ","), **tokens_pt = tokens;
      while(tok != NULL) {
        switch(format[pos]) {
          case 's':
            temp.push(tok);
            break;
          case 'd':
            d = strtod(tok, NULL);
            temp.push(&d);
            break;
          case 'x':
            if( strcmp(tok, *tokens_pt++) )
              goto SKIPLINE;
            break;
          case '\0':
            goto SKIPLINE;
        }
        tok = strtok(NULL, ",");
        pos++;
      }

      if(pos < d_len) {
        goto SKIPLINE;
      }
      data_vect.push_back(temp);
      return;

SKIPLINE:
      temp.free_mem();
      return;
    }
};

int main (int argc, char *argv[]) {
  clock_t start_time = clock();

  if(argc < 3) {
    std::cerr << "Expected at least 2 arguments: csvParse FILENAME FORMAT tokens..." << std::endl;
    return 1;
  }

  char *fn = argv[1];
  char *check = argv[2];
  int count = 0;
  do {
    if( *check != 'x' && *check != 'd' && *check != 's') {
      std::cerr << "Illegal FORMAT identifier. Please use only 'd' double 's' string and 'x' token." << std::endl;
      return 1;
    }
    if(*check == 'x')
      count++;
  } while(*++check != '\0');
  if(count != (argc - ARGS_OFFSET)) {
    std::cerr << "The number of tokens passed as arguments does not match the identifiers in FORMAT" << std::endl;
    return 1;
  }
  char **tokens = new char*[count];
  for(int i = 0; i < count; i++)
    tokens[i] = argv[i + ARGS_OFFSET];

  CSV::initialize(fn, argv[2], tokens);

  CSV::parser temp;
  int line_count = 0;
  while( !temp.readline() ){ line_count++; }

  std::cout << "Finished reading CSV file. Took " << (double) (clock() - start_time) / CLOCKS_PER_SEC * 1000 << " milliseconds." << std::endl;
  std::cout << "Found " << temp.data_vect.size() << " matches in " << line_count << " lines." << std::endl;
  std::cout << "Consolidating CSV::data vector columns..." << std::endl;
  //char **test = (char **) temp.consolidate();
  //for(int i = 0; i < 3; i++)
  //  std::cout << test[i] << std::endl;

  temp.consolidate();
  temp.consolidate();

  std::cout << "Consolidation complete. Total runtime was " << (double) (clock() - start_time) / CLOCKS_PER_SEC * 1000 << " milliseconds. Have a nice day!" << std::endl;

  CSV::exit();

  return 0;
}
