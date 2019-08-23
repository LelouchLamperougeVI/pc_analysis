#ifndef __CSV__H
#define __CSV__H

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

namespace CSV {
  extern std::ifstream stream; //file stream
  extern char *format; //format identifier
  extern char **tokens; //search tokens
  extern std::mutex mute; //mutex guard
  extern size_t d_size; //the size of each datum on memory
  extern int d_len; //number of elements stored in each datum

  void initialize(std::string &fn, char *iformat, char **itokens);
  void exit();

  class parser;
  class data;
};

class CSV::data {
  public:
    data();
    void memfree();
    void push(void *new_data);
    void *pop();

  private:
    void *datum;
    void *write_head;
    char *currPos;
};

class CSV::parser {
  public:
    parser();
    ~parser();
    int size();
    int readline();
    void *consolidate(char &flag);

  private:
    std::vector<data> data_vect;
    std::vector<char **> char_col_vect;
    std::vector<double *> double_col_vect;
    char buff[BUFF_LEN];
    char *sp = format; //old habits from assembly will never change...

    void parsel(char *line);
};
#endif
