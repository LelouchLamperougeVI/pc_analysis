#include "MATW.h"

MATW::writer::writer(std::string &fn, char *of, mxClassID type, int *size) {
  const char *fn_char = fn.c_str();
  matfile = matOpen(fn_char, "w");
  if(matfile == NULL)
    throw "Failed to open .mat file for writing.";

  switch(type) {
    case mxDOUBLE_CLASS: {
      dst_type = 'd';
      dst = mxCreateDoubleMatrix(size[0], size[1], mxREAL);
      break;
              }
    case mxCELL_CLASS: {
      dst_type = 'c';
      dst = mxCreateCellMatrix(size[0], size[1]);
      break;
              }
    default:
        throw "Must specify type of destination mxArray as either 'd' for double or 'c' for cell.";
  };

  out_fn = of;
  M = mxGetM(dst);
  N = mxGetN(dst);
}

MATW::writer::~writer() {
  matPutVariable(matfile, out_fn, dst);
  matClose(matfile);
  for(int i = 0; i < buff.size(); i++)
    for(int j = 0; j < M; j++)
      mxDestroyArray(*buff[i]++);
  buff.clear();
}

void MATW::writer::write(void *src, char type, int index) {
  if( dst_type == 'd' && type == 'd' )
    double2double(src, index);
  else if( dst_type == 'c' && type == 'd' )
    double2cell(src, index);
  else if( dst_type == 'c' && type == 's' )
    char2cell(src, index);
  else
    throw "Bad source.";
}

void MATW::writer::char2cell(void *src, int col) {
  char **pt = (char **) src;
  mxArray **temp = new mxArray*[M];
  for(int i = 0; i < M; i++) {
    temp[i] = mxCreateString(*pt++);
    mxSetCell(dst, i + col*M, temp[i]);
  }
  buff.push_back(temp);
}

void MATW::writer::double2cell(void *src, int col) {
  double *pt = (double *) src;
  mxArray **temp = new mxArray*[M];
  for(int i = 0; i < M; i++) {
    temp[i] = mxCreateDoubleScalar(*pt++);
    mxSetCell(dst, i + col*M, temp[i]);
  }
  buff.push_back(temp);
}

void MATW::writer::double2double(void *src, int col) {
  memcpy( (void *) (mxGetPr(dst) + col*M), src, M * sizeof(double) );
}
