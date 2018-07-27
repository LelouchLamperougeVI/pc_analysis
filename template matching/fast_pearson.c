/** Faster pearson correlation
   rho = fast_pearson(x, y)
   Y is a single column vector
 */

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "pthread.h"

#define NUMTHREADS 7

struct ccData {
        int m, idx;
        int *numTasks;
        double *mu_x, *mu_y, *X, *Y, *cc;
};

void getMean(double *X, double *mean, int m){
        int i;
        for(i = 0; i < m; i++)
              *mean += X[i];
        *mean /= m;
}

void subMean(double *X, double *mean, int m){
        int i;
        for(i = 0; i < m; i++)
              X[i] -= *mean;
}

void corrcoef(struct ccData *data){
        int m = data->m;
        int *numTasks = data->numTasks;
        int idx = data->idx;
        double *mu_x = data->mu_x;
        double *mu_y = data->mu_y;
        double *X = data->X;
        double *Y = data->Y;
        double *cc = data->cc;

        double top, left, right;

        int i, j;
        int block_i = 0;
        for(i = 0; i < idx; i++)
                block_i += numTasks[i];

        for(i = 0; i < numTasks[idx]; i++){
                top = 0.0; left = 0.0; right = 0.0;
                getMean(X + ((block_i+i)*m), mu_x + block_i+i, m);
                for(j = 0; j < m; j++){
                        top += (X[(block_i+i)*m + j] - mu_x[block_i+i]) * Y[j];
                        left += (X[(block_i+i)*m + j] - mu_x[block_i+i]) * (X[(block_i+i)*m + j] - mu_x[block_i+i]);
                        right += Y[j] * Y[j];
                }
                cc[block_i+i] = top / sqrt(left * right);
        }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        if(nrhs > 2)
                mexErrMsgTxt("Too many inputs");
        if(nrhs < 2)
                mexErrMsgTxt("Missing inputs");
        if(nlhs > 1)
                mexErrMsgTxt("Too many outputs assigned");

        double *mu_x, *mu_y, *X, *Y, *cc;
        mxArray *mu_x_pr, *mu_y_pr;
        int i, m, n;
        int tasksPerThread[NUMTHREADS] = { 0 };

        X = mxGetPr(prhs[0]);
        Y = mxGetPr(prhs[1]);
        m = mxGetM(prhs[0]);
        n = mxGetN(prhs[0]);

        plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
        cc = mxGetPr(plhs[0]);
        mu_x_pr = mxCreateDoubleMatrix(n, 1, mxREAL);
        mu_x = mxGetPr(mu_x_pr);
        mu_y_pr = mxCreateDoubleScalar(0.0);
        mu_y = mxGetPr(mu_y_pr);

        getMean(Y, mu_y, m);
        subMean(Y, mu_y, m);

        struct ccData data[NUMTHREADS];
        pthread_t thread[NUMTHREADS];

        for(i = 0; i < n; i++)
                tasksPerThread[i % NUMTHREADS]++;

        for(i = 0; i < NUMTHREADS; i++) {
                data[i].X = X;
                data[i].Y = Y;
                data[i].mu_x = mu_x;
                data[i].mu_y = mu_y;
                data[i].cc = cc;
                data[i].m = m;
                data[i].numTasks = tasksPerThread;
                data[i].idx = i;
        }

        for(i = 0; i < NUMTHREADS; i++)
                pthread_create(&thread[i], NULL, corrcoef, &data[i]);

        for(i = 0; i < NUMTHREADS; i++)
                pthread_join(thread[i], NULL);

        mxDestroyArray(mu_x_pr);
        mxDestroyArray(mu_y_pr);
        return;
}
