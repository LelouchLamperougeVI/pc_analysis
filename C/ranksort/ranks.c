#include "ranks.h"

void *getRanks(void *arg){
        struct ranksData *data = (struct ranksData *) arg;
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
                        while(A[(block_i + i + 1)*m - buff - 1] == A[(block_i + i + 1)*m - buff - 2] && buff < (m - 1)) {
                                buff++;
                        }
                        do {
                                ranks[(block_i + i + 1)*m - k - 1] = (buff - count) / 2.0 + count + 1.0;
                                k++;
                        } while(k <= buff);
                        count = buff + 1.0;
                } while(count < m);
        }
}
