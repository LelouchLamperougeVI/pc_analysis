#ifndef __MATW__
#define __MATW__

#include <vector>
#include <cstring>
#include "mat.h"

namespace MATW { //mat file writer
  class writer;
};

class MATW::writer {
  public:
    writer(char *fn, char *of, mxClassID type, int *size);
    ~writer();
    void write(void *src, char type, int index);

  private:
    MATFile *matfile;
    char dst_type, src_type, *out_fn;
    mwSize M, N;
    mxArray *dst;
    std::vector<mxArray **> buff;

    void char2cell(void *src, int col);
    void double2cell(void *src, int col);
    void double2double(void *src, int col);
};


#endif
