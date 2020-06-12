#include "CSV.h"

std::ifstream CSV::stream;
char *CSV::format;
char **CSV::tokens;
std::mutex CSV::mute;
size_t CSV::d_size;
int CSV::d_len;

void CSV::initialize(std::string &fn, char *iformat, char **itokens) {
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

/********************CSV::data********************/
CSV::data::data() {
  datum = malloc(d_size);
  write_head = datum;
  currPos = format;
}

void CSV::data::memfree() { //lol most pathetic memory leak ever...
  free(datum);
}

void CSV::data::push(void *new_data) {
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

void *CSV::data::pop() {
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

/********************CSV::parser********************/
CSV::parser::parser() {
  while(*++sp != '\0') {}
}

CSV::parser::~parser() {
  for(int i = 0; i < data_vect.size(); i++)
    data_vect[i].memfree();
  data_vect.clear();

  for(int i = 0; i < char_col_vect.size(); i++)
    delete [] char_col_vect[i];
  for(int i = 0; i < double_col_vect.size(); i++)
    delete [] double_col_vect[i];
  char_col_vect.clear();
  double_col_vect.clear();
}

int CSV::parser::size() {
  return data_vect.size();
}

int CSV::parser::readline() {
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

void *CSV::parser::consolidate(char &flag) { //pop out column from last data element
    SKIP: switch(*--sp) {
    case 'd': {
      double *column = new double[data_vect.size()];
      for(int i = 0; i < data_vect.size(); i++) {
        double *temp = (double *) data_vect[i].pop();
        column[i] = *temp;
      }
      double_col_vect.push_back(column);
      flag = *sp;
      return column;
    }
    case 's': {
      char **column = new char*[data_vect.size()];
      for(int i = 0; i < data_vect.size(); i++) {
        column[i] = new char[CHAR_LEN];
        column[i] = (char *) data_vect[i].pop();
      }
      char_col_vect.push_back(column);
      flag = *sp;
      return column;
    }
    case 'x':
      goto SKIP;
    default:
      return NULL;
  }
}

void CSV::parser::parsel(char *line) { //parse line
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
  temp.memfree();
  return;
}
