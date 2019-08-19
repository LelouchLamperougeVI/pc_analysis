/*************************************************************************
 * MATLAB MEX ROUTINE csv_parse.c
 *
 * Index keywords in Smoothwalk CSV file
 * 1: userEntry
 * 2: position
 * 3: count-21
 * 4: pickup
 *
 * Direct complaints to HaoRan Chang
 ************************************************************************/

#include "string.h"
#include "mex.h"
#include "pthread.h"

#define NUMTHREADS 8

struct threadData{
    int start, stop;
    const mxArray *cellPr;
    double *idx;
};

void *findIdx(struct threadData *data){
    int start = data->start;
    int stop = data->stop;
    const mxArray *cellPr = data->cellPr;
    double *idx = data->idx;

    const mxArray *str;

    char *buff;
    int buflen = 100;

    buff = mxMalloc(buflen);

    int i;
    for(i=start; i<stop; i++){
        str = mxGetCell(cellPr,i);
        mxGetString(str, buff, buflen);
        if(!strcmp("parse", buff)){ //start
            idx[i]=1;
        }else if(!strcmp("pin-state", buff)){
            idx[i]=2;
        }else if(!strcmp("sprite", buff)){ //black
            idx[i]=3;
        }else if(!strcmp("object", buff)){ //object
            idx[i]=4;
        }else{
            idx[i]=0;
        }
    }

    mxFree(buff);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    int numElements, i;

    const mxArray *cellPr = prhs[0];
    numElements = mxGetNumberOfElements(cellPr);

    plhs[0] = mxCreateDoubleMatrix(1, numElements, mxREAL);
    double *idx = mxGetPr(plhs[0]);

    struct threadData data[NUMTHREADS];
    pthread_t thread[NUMTHREADS];
    int tasksPerThread = numElements / NUMTHREADS;

    for(i=0; i<NUMTHREADS; i++){
        data[i].start = i*tasksPerThread;
        data[i].stop = (i+1)*tasksPerThread;
        data[i].cellPr = cellPr;
        data[i].idx = idx;
    }
    data[NUMTHREADS - 1].stop = numElements;

    for(i=0; i<NUMTHREADS; i++){
        pthread_create(&thread[i], NULL, findIdx, &data[i]);
    }

    for(i=0; i<NUMTHREADS; i++){
        pthread_join(thread[i], NULL);
    }

    return;
}
