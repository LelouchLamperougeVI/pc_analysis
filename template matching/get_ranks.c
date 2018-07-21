/** Descending fractional ranking
 */

#include "mex.h"
#include "matrix.h"
#include "pthread.h"

#define NUMTHREADS 7

struct threadData {
        int m, idx;
        int *numTasks;
        double *A, *ranks;
};

void getRanks(struct threadData *data){
        int m = data->m;
        int *numTasks = data->numTasks;
        int idx = data->idx;
        double *A = data->A;
        double *ranks = data->ranks;

        int i, block_i, buff, k;
        double count;

        block_i = 0;
        for(i = 0; i < idx; i++)
                block_i += numTasks[i];

        for(i = 0; i < numTasks[idx]; i++) {
                count = 0.0;
                do {
                        k = (int) count;
                        buff = (int) count;
                        while(A[(block_i + i)*m - buff - 1] == A[(block_i + i)*m - buff - 2] && buff < (m - 1))
                                buff++;
                        do {
                                ranks[(block_i + i)*m - k - 1] = (buff - count) / 2.0 + count + 1.0;
                                k++;
                        } while(k <= buff);
                        count = buff + 1.0;
                } while(count < m);
        }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        if(nrhs > 1)
                mexErrMsgTxt("Too many inputs");
        if(nrhs < 1)
                mexErrMsgTxt("No input given");
        if(nlhs > 1)
                mexErrMsgTxt("Too many outputs assigned");

        double *A, *ranks;
        int i, m, n;
        int tasksPerThread[NUMTHREADS] = { 0 };

        A = mxGetPr(prhs[0]);
        m = mxGetM(prhs[0]);
        n = mxGetN(prhs[0]);

        plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
        ranks = mxGetPr(plhs[0]);

        struct threadData data[NUMTHREADS];
        pthread_t thread[NUMTHREADS];

        for(i = 0; i < n; i++)
                tasksPerThread[i % NUMTHREADS]++;

        for(i = 0; i < NUMTHREADS; i++)
                mexPrintf("%d ", tasksPerThread[i]);

        for(i = 0; i < NUMTHREADS; i++) {
                data[i].A = A;
                data[i].ranks = ranks;
                data[i].m = m;
                data[i].numTasks = tasksPerThread;
                data[i].idx = i;
        }

        for(i = 0; i < NUMTHREADS; i++)
                pthread_create(&thread[i], NULL, getRanks, &data[i]);

        for(i = 0; i < NUMTHREADS; i++)
                pthread_join(thread[i], NULL);

        return;
}
