/** Helper function for the k-NN estimator of mutual information (Kraskov et al. 2004)
   n = Iknn(A,k)
 */

#include "mex.h"
#include "matrix.h"
#include "pthread.h"
#include "ranks.h"
#include "sort.h"
#include "string.h"

#define NUMTHREADS 7

struct knnData {
        int m, t_idx;
        int *numTasks;
        double *idx, *ranks;
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        if(nrhs > 2)
                mexErrMsgTxt("Too many inputs");
        if(nrhs < 1)
                mexErrMsgTxt("No input given");
        if(nlhs > 2)
                mexErrMsgTxt("Too many outputs assigned");

        double *A, *sorted, *idx, *ranks;
        int i, m, n;
        int tasksPerThread[NUMTHREADS] = { 0 };

        A = mxGetPr(prhs[0]);
        m = mxGetM(prhs[0]);
        n = mxGetN(prhs[0]);

        plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
        sorted = mxGetPr(plhs[0]);
        plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
        idx = mxGetPr(plhs[1]);
        plhs[2] = mxCreateDoubleMatrix(m, n, mxREAL);
        ranks = mxGetPr(plhs[2]);


        for(i = 0; i < m; i++){
          A[i]+=
        }


        struct sortData data[NUMTHREADS];
        pthread_t thread[NUMTHREADS];

        for(i = 0; i < n; i++)
                tasksPerThread[i % NUMTHREADS]++;

        for(i = 0; i < NUMTHREADS; i++) {
                data[i].A = A;
                data[i].sorted = sorted;
                data[i].m = m;
                data[i].numTasks = tasksPerThread;
                data[i].idx = idx;
                data[i].t_idx = i;
        }

        for(i = 0; i < NUMTHREADS; i++)
                pthread_create(&thread[i], NULL, quicksort, &data[i]);

        for(i = 0; i < NUMTHREADS; i++)
                pthread_join(thread[i], NULL);


        //ranks
        struct ranksData rData[NUMTHREADS];
        pthread_t rThread[NUMTHREADS];

        for(i = 0; i < NUMTHREADS; i++) {
                rData[i].A = sorted;
                rData[i].ranks = ranks;
                rData[i].m = m;
                rData[i].numTasks = tasksPerThread;
                rData[i].idx = i;
        }

        for(i = 0; i < NUMTHREADS; i++)
                pthread_create(&rThread[i], NULL, getRanks, &rData[i]);

        for(i = 0; i < NUMTHREADS; i++)
                pthread_join(rThread[i], NULL);

        //sort ranks
        struct knnData iData[NUMTHREADS];
        pthread_t iThread[NUMTHREADS];

        for(i = 0; i < NUMTHREADS; i++) {
                iData[i].idx = idx;
                iData[i].ranks = ranks;
                iData[i].m = m;
                iData[i].numTasks = tasksPerThread;
                iData[i].t_idx = i;
        }

        for(i = 0; i < NUMTHREADS; i++)
                pthread_create(&iThread[i], NULL, indexRanks, &iData[i]);

        for(i = 0; i < NUMTHREADS; i++)
                pthread_join(iThread[i], NULL);

        // indexRanks(ranks, idx, m, n);

        return;
}
