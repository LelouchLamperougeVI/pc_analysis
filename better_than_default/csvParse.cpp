/* csvParse.cpp
 * Takes a csv file, parse through it with the specified format identifier
 * and save the results as a .mat file for use in MATLAB
 *
 * Compilation: mex -client engine csvParse.cpp CSV.cpp MATW.cpp
 * Compiled with g++ 9.1.0 on MATLAB R2019a and tested on Arch Linux kernel 5.2.8
 *
 * Usage: csvParse [FILENAME] [FORMAT] [tokens...]
 *  FILENAME    - Input CSV file, given with full path if not in pwd.
 *  FORMAT      - The template of the lines to be parsed. Available identifiers
 *                are 'd', 's' and 'x'.
 *  [tokens]    - For each 'x' identifier, provide a search token.
 *
 * Output: a .mat file with the same name as the input file.
 *
 * e.g. In a typical two-photon imaging session, one might want to extract all the
 * entries that are formatted as
 *     'time-stamp,pin number,XX,pin-state,XX,object,"elephant"'.
 * In this case, the format identifiers should be given as 'dxdxdxs', which inform
 * the program to extract three 'doubles' and one string. 'x's dictate the locations
 * of search tokens. The location of each identifier is relative to tokens on a
 * line separated by comma delimiters. The full command is
 *     'csvParse example.csv dxdxdxs pin\ number pin-state object'.
 *
 * TODO: code is bloated... but it works... not the most elegant solution...
 * TODO: Originally was gonna use threads, but the speed was so much faster than what
 *       MATLAB could achieve just off a single thread. So no need for that
 *       functionality anymore. Remove.
 *
 * By HaoRan Chang, Ph.D. candidate
 * Canadian Centre for Behavioural Neuroscience,
 * Department of Neuroscience, University of Lethbridge,
 * Lethbridge, Alberta, Canada
 */

#include "MATW.h"
#include "CSV.h"
#include <iostream>
#include <ctime>

#define ARGS_OFFSET 3

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

  CSV::parser parser;
  int line_count = 0;
  while( !parser.readline() ){ line_count++; }

  std::cout << "Finished reading CSV file. Took " << (double) (clock() - start_time) / CLOCKS_PER_SEC * 1000 << " milliseconds." << std::endl;
  std::cout << "Found " << parser.size() << " matches in " << line_count << " lines." << std::endl;
  std::cout << "Consolidating CSV::data vector columns..." << std::endl;

  char *tok = strtok(fn, "/\\");
  char *of;
  char extension[] = ".mat";
  while(tok != NULL) {
    of = tok;
    tok = strtok(NULL, "/\\");
  }
  tok = strtok(of, ".");
  of = new char[CHAR_LEN];
  strcpy(of, tok);
  strcpy(of + strlen(tok), extension);

  int dims[2] = {parser.size(), CSV::d_len};
  char *col_pt = CSV::format;
  char varname[] = "CSV";
  MATW::writer *writer;
  while(*col_pt != '\0' && *col_pt != 's') {col_pt++;}
  switch(*col_pt) {
    case 's': {
      writer = new MATW::writer(of, varname, mxCELL_CLASS, dims);
      break;
    }
    default:
      writer = new MATW::writer(of, varname, mxDOUBLE_CLASS, dims);
      break;
  }

  char flag;
  void *ret;
  count = CSV::d_len-1;
  do {
    ret = parser.consolidate(flag);
    switch(flag) {
      case 'd': {
        writer->write((double *) ret, flag, count--);
        break;
      }
      case 's': {
        writer->write((char **) ret, flag, count--);
        break;
      }
    }
  } while(count > -1);

  std::cout << "Consolidation complete. Total runtime was " << (double) (clock() - start_time) / CLOCKS_PER_SEC * 1000 << " milliseconds. Have a nice day!" << std::endl;

  delete [] tokens;
  CSV::exit();
  delete writer;
  delete [] of;

  return 0;
}
