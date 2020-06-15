#include <stdio.h>
#include <stdlib.h>
#include "propagation-sphere.h"

// Define global variables here
int TIME_ORDER = 2;
int DIMS = 2;

int main(int argc, char const *argv[])
{

    int BORDER_SIZE = 0;
    int SPACE_ORDER = 2;
    int time_m = 1;
    int time_M = 8;
    int GRID_SIZE = 44;
    int x_m = (int)BORDER_SIZE + SPACE_ORDER;
    int x_M = (int)BORDER_SIZE + SPACE_ORDER + GRID_SIZE;
    int y_m = (int)BORDER_SIZE + SPACE_ORDER;
    int y_M = (int)BORDER_SIZE + SPACE_ORDER + GRID_SIZE;

    int size_u[] = {GRID_SIZE + 2 * BORDER_SIZE + 2 * SPACE_ORDER, GRID_SIZE + 2 * BORDER_SIZE + 2 * SPACE_ORDER};

    float **vp;
    float ***u;

    // Allocate data
    printf("Allocating u array with dimensions %dx%dx%d...\n", TIME_ORDER + 1, size_u[0], size_u[1]);
    u = (float ***)malloc(sizeof(float **) * (TIME_ORDER + 1));
    for (int i = 0; i < TIME_ORDER + 1; i++)
    {
        u[i] = (float **)malloc(sizeof(float *) * size_u[0]);
        for (int j = 0; j < size_u[0]; j++)
        {
            u[i][j] = (float *)malloc(sizeof(float) * size_u[1]);
        }
    }

    printf("Allocating vp array with dimensions %dx%d...\n", size_u[0], size_u[1]);
    vp = (float **)malloc(sizeof(float *) * size_u[0]);
    for (int j = 0; j < size_u[0]; j++)
    {
        vp[j] = (float *)malloc(sizeof(float) * size_u[1]);
    }

    // Initializing values
    printf("Initializing data...\n");
    for (int j = 0; j < size_u[0]; j++)
    {
        for (int k = 0; k < size_u[1]; k++)
        {
            u[0][j][k] = 0.0;
            u[1][j][k] = 0.0;
            u[2][j][k] = 0.0;
            vp[j][k] = 1.5;
        }
    }
    // Source injection
    printf("Source injection...\n");
    u[0][22][22] = 1.;

    // Print intial data
    print_array_2d(u[0], size_u[0], size_u[1]);

    printf("Wave propagation...\n");
    float r1 = 0.0784;
    for (int time = time_m, t0 = (time) % (3), t1 = (time + 1) % (3), t2 = (time + 2) % (3);
         time <= time_M;
         time += 1, t0 = (time) % (3), t1 = (time + 1) % (3), t2 = (time + 2) % (3))
    {
        printf("t0=%d,t1=%d,t2=%d", t0, t1, t2);
        for (int x = x_m / 2; x < x_M; x += 1)
        {
            for (int y = y_m; y < y_M; y += 1)
            {
                //float r0 = vp[(x*(size_u[0] + 2)) + y + 2] * vp[(x*(size_u[0] + 2)) + y + 2];
                printf("[%d][%d] ", x + 3, y + 3);
                float r0 = vp[x + 2][y + 2] * vp[x + 2][y + 2];
                u[t1][x + 2][y + 2] = -4.0F * r0 * r1 * u[t0][x + 2][y + 2] +
                                      1.0F * (r0 * r1 * u[t0][x + 1][y + 2] +
                                              r0 * r1 * u[t0][x + 2][y + 1] +
                                              r0 * r1 * u[t0][x + 2][y + 3] +
                                              r0 * r1 * u[t0][x + 3][y + 2] -
                                              u[t2][x + 2][y + 2]) +
                                      2.0F * u[t0][x + 2][y + 2];
            }
            printf("\n");
        }
    }

    // Print result
    print_array_2d(u[1], size_u[0], size_u[1]);

    return 0;
}

void print_array_2d(float **u, int x_size, int y_size)
{

    for (int j = 0; j < x_size; j++)
    {
        for (int k = 0; k < y_size; k++)
        {
            printf("%.7f ", u[k][j]);
        }
        printf("\n");
    }
    printf("\n\n");
}
