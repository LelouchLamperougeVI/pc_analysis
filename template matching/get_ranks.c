/** Fractional ranking algortihm
   ranks = get_ranks(sorted_matrix)
 */

#include "mex.h"
#include "matrix.h"
#include "pthread.h"
#include "ranks.h"

#define NUMTHREADS 7

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

        struct ranksData data[NUMTHREADS];
        pthread_t thread[NUMTHREADS];

        for(i = 0; i < n; i++)
                tasksPerThread[i % NUMTHREADS]++;

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
