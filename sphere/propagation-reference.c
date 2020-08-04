#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "propagation-sphere.h"
#include "source.h"

// Define global variables here
int TIME_ORDER = 2;

int main(int argc, char const *argv[])
{

    int BORDER_SIZE = 0;
    int SPACE_ORDER = 2;
    int time_m = 1;
    int time_M = 900;
    int GRID_SIZE = 10000;
    int x_m = (int)BORDER_SIZE + SPACE_ORDER;
    int x_M = (int)BORDER_SIZE + SPACE_ORDER + GRID_SIZE;
    int y_m = (int)BORDER_SIZE + SPACE_ORDER;
    int y_M = (int)BORDER_SIZE + SPACE_ORDER + GRID_SIZE;

    int size_u[] = {GRID_SIZE + 2 * BORDER_SIZE + 2 * SPACE_ORDER,
                    GRID_SIZE + 2 * BORDER_SIZE + 2 * SPACE_ORDER};

    float **vp;
    float ***u;

    // Allocate data
    // printf("Allocating u array with dimensions %d x %d x %d...\n", TIME_ORDER + 1, size_u[0], size_u[1]);
    u = (float ***)malloc(sizeof(float **) * (TIME_ORDER + 1));
    for (int i = 0; i < TIME_ORDER + 1; i++)
    {
        u[i] = (float **)malloc(sizeof(float *) * size_u[0]);
        for (int j = 0; j < size_u[0]; j++)
        {
            u[i][j] = (float *)malloc(sizeof(float) * size_u[1]);
        }
    }

    // printf("Allocating vp array with dimensions %d x %d...\n", size_u[0], size_u[1]);
    vp = (float **)malloc(sizeof(float *) * size_u[0]);
    for (int j = 0; j < size_u[0]; j++)
    {
        vp[j] = (float *)malloc(sizeof(float) * size_u[1]);
    }

    // Initializing values
    // printf("Initializing data...\n");
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
    // Ricker source
    float fpeak = 1;
    float *source = ricker_source(fpeak, time_M - time_m + 1, 1); // dt é o incremento no laço do tempo...
    int source_location[2] = {size_u[0] / 2, size_u[1] / 2};      // Middle_point of each dimension.

    // Print intial data
    // print_array_2d(u[0], size_u[0], size_u[1]);

    // printf("Fonte: \n");
    // for (int i = 0; i < time_M - time_m + 1; i++)
    // {
    //     printf("%.3f ", source[i]);
    // }
    // printf("\n");

    // printf("Wave propagation...\n");
    float r1 = 0.0784;
    for (int time = time_m, t0 = (time) % (3), t1 = (time + 1) % (3), t2 = (time + 2) % (3);
         time <= time_M;
         time += 1, t0 = (time) % (3), t1 = (time + 1) % (3), t2 = (time + 2) % (3))
    {

        for (int x = x_m - 1; x < x_M - 1; x += 1)
        {
            for (int y = y_m - 1; y < y_M - 1; y += 1)
            {
                float r0 = vp[x + 2][y + 2] * vp[x + 2][y + 2];
                u[t1][x + 2][y + 2] = -4.0F * r0 * r1 * u[t0][x + 2][y + 2] +
                                      1.0F * (r0 * r1 * u[t0][x + 1][y + 2] +
                                              r0 * r1 * u[t0][x + 2][y + 1] +
                                              r0 * r1 * u[t0][x + 2][y + 3] +
                                              r0 * r1 * u[t0][x + 3][y + 2] -
                                              u[t2][x + 2][y + 2]) +
                                      2.0F * u[t0][x + 2][y + 2];
            }
        }

        // Inject source
        u[t1][source_location[0]][source_location[1]] += source[time - 1];
    }

    // Print result
    // print_array_2d(u[1], size_u[0], size_u[1]);

   // Free resources
    free(source);
    free(vp);
    free(u[0]);
    free(u[1]);
    free(u[2]);

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