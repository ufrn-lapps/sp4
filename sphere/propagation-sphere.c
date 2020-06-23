#include <stdio.h>
#include <stdlib.h>
#include "propagation-sphere.h"

#define max(a, b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a, b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
// Define global variables here
int TIME_ORDER = 2;

int main(int argc, char const *argv[])
{

    int BORDER_SIZE = 0;
    int SPACE_ORDER = 2;
    int time_m = 1;
    int time_M = 8;

    int grid_points[2] = {44, 44};   // Number of points in each dimension of the grid.
    int points_distance[2] = {1, 1}; // Distance between each point of the grid in each dimension. in meters.
                                     // The grid length will be given by grid_points * points_distance

    int x_m = (int)BORDER_SIZE + SPACE_ORDER;
    int x_M = (int)BORDER_SIZE + SPACE_ORDER + grid_points[0];
    int y_m = (int)BORDER_SIZE + SPACE_ORDER;
    int y_M = (int)BORDER_SIZE + SPACE_ORDER + grid_points[1];

    int size_u[] = {grid_points[0] + 2 * BORDER_SIZE + 2 * SPACE_ORDER,  // The array u, is composed by the grid points + 
                    grid_points[1] + 2 * BORDER_SIZE + 2 * SPACE_ORDER}; // border on the left and right + space order cells left and right.

    float **vp;
    float ***u;

    // Allocate data
    printf("Allocating u array with dimensions %d x %d x %d...\n", TIME_ORDER + 1, size_u[0], size_u[1]);
    u = (float ***)malloc(sizeof(float **) * (TIME_ORDER + 1));
    for (int i = 0; i < TIME_ORDER + 1; i++)
    {
        u[i] = (float **)malloc(sizeof(float *) * size_u[0]);
        for (int j = 0; j < size_u[0]; j++)
        {
            u[i][j] = (float *)malloc(sizeof(float) * size_u[1]);
        }
    }

    printf("Allocating vp array with dimensions %d x %d...\n", size_u[0], size_u[1]);
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
    // Max velocity
    float vp_max = 1.5; // In this case, we only have this velocity, but if we have many, we get the greater.

    // Source injection
    printf("Source injection...\n");
    int source_location[2] = {size_u[0] / 2, size_u[1] / 2}; // Middle_point of each dimension.
    // Injecting 1 in the source location.
    printf("\tSource Location: %d %d\n", source_location[0], source_location[1]);
    u[0][source_location[0]][source_location[1]] = 1.;

    // Print intial data
    // print_array_2d(u[0], size_u[0], size_u[1]);

    printf("Wave propagation...\n");
    //printf("\tPropagation time: %d until %d\n", time_m, time_M);
    //printf("\tWave front velocity: %.2f\n", vp_max);
    float r1 = 0.0784;
    int x_init, x_end;
    int y_init, y_end;
    for (int time = time_m, t0 = (time) % (3), t1 = (time + 1) % (3), t2 = (time + 2) % (3);
         time <= time_M;
         time += 1, t0 = (time) % (3), t1 = (time + 1) % (3), t2 = (time + 2) % (3))
    {

        x_init = max((x_m - 1), (int)(source_location[0] - ((vp_max * time) / points_distance[0])));
        x_end = min((x_M - 1), (int)(source_location[0] + ((vp_max * time) / points_distance[0])));
        // printf("\tX range: from: %d until %d\n", x_init, x_end);

        y_init = max((y_m - 1), (int)(source_location[1] - ((vp_max * time) / points_distance[1])));
        y_end = min((y_M - 1), (int)(source_location[1] + ((vp_max * time) / points_distance[1])));
        //printf("\tY range: from: %d until %d\n", y_init, y_end);

        for (int x = x_init; x < x_end; x += 1)
        {
            for (int y = y_init; y < y_end; y += 1)
            {
                //float r0 = vp[(x*(size_u[0] + 2)) + y + 2] * vp[(x*(size_u[0] + 2)) + y + 2];
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
    }

    // Print result
    // print_array_2d(u[1], size_u[0], size_u[1]);

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
