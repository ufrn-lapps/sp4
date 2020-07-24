#include <stdio.h>
#include <stdlib.h>
#include "math.h"
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

    int grid_points[2] = {44, 44};         // Number of points in each dimension of the grid.
    float points_distance[2] = {1.0, 1.0}; // Distance between each point of the grid in each dimension. in meters.
                                           // The grid length will be given by grid_points * points_distance

    int x_m = (int)BORDER_SIZE + SPACE_ORDER;
    int x_M = (int)BORDER_SIZE + SPACE_ORDER + grid_points[0];
    int y_m = (int)BORDER_SIZE + SPACE_ORDER;
    int y_M = (int)BORDER_SIZE + SPACE_ORDER + grid_points[1];

    int size_u[] = {grid_points[0] + 2 * BORDER_SIZE + 2 * SPACE_ORDER,  // The array u, is composed by the grid points +
                    grid_points[1] + 2 * BORDER_SIZE + 2 * SPACE_ORDER}; // border on the left and right + space order cells left and right.

    float **vp;
    float ***u;
    // float u[3][48][48];

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
    // Max velocity
    float vp_max = 1.5; // In this case, we only have this velocity, but if we have many, we must get the greater.

    // Source injection
    // printf("Source injection...\n");
    int source_location[2] = {size_u[0] / 2, size_u[1] / 2}; // Middle_point of each dimension.
    // Injecting 1 in the source location.
    // printf("\tSource Location: %d %d\n", source_location[0], source_location[1]);
    u[time_m % 3][source_location[0]][source_location[1]] = 1.;
    // printf("(%d,%d)=%f ", source_location[0], source_location[1], u[0][source_location[0]][source_location[1]]);

    // Print intial data
    // print_array_2d(u[0], size_u[0], size_u[1]);

    // printf("Wave propagation...\n");
    //printf("\tPropagation time: %d until %d\n", time_m, time_M);
    //printf("\tWave front velocity: %.2f\n", vp_max);
    float r1 = 0.0784;
    float x_min, x_max;
    float y_min, y_max;
    float raio = 0.0;

    for (int time = time_m, t0 = (time) % (3), t1 = (time + 1) % (3), t2 = (time + 2) % (3);
         time <= time_M;
         time += 1, t0 = (time) % (3), t1 = (time + 1) % (3), t2 = (time + 2) % (3))
    {
        // Calculate dimension extremes
        raio = vp_max * time; // Em metros

        y_min = max((y_m - 1) * points_distance[1], source_location[1] * points_distance[1] - (raio)); // Em metros

        y_max = min((y_M - 1) * points_distance[1], source_location[1] * points_distance[1] + (raio)); // Em metros

        for (int y = (int)floor(y_min / points_distance[1]); y < (int)ceil(y_max / points_distance[1]); y += 1) // Converte para pontos cartesianos.
        {
            // Em metros
            x_min = max((x_m - 1) * points_distance[0], source_location[0] * points_distance[0] - sqrt(raio * raio - pow((y - source_location[1] * points_distance[1]), 2)));
            // Em metros
            x_max = min((x_M - 1) * points_distance[0], source_location[0] * points_distance[0] + sqrt(raio * raio - pow((y - source_location[1] * points_distance[1]), 2)));

            for (int x = (int)floor(x_min / points_distance[0]); x < (int)ceil(x_max / points_distance[0]); x += 1) // Converte para pontos cartesianos.
            {

                float r0 = vp[x][y] * vp[x][y];
                u[t1][x][y] = -4.0F * r0 * r1 * u[t0][x][y] +
                              1.0F * (r0 * r1 * u[t0][x - 1][y] +
                                      r0 * r1 * u[t0][x][y - 1] +
                                      r0 * r1 * u[t0][x][y + 1] +
                                      r0 * r1 * u[t0][x + 1][y] -
                                      u[t2][x][y]) +
                              2.0F * u[t0][x][y];
            }
        }
        print_array_2d(u[1], size_u[0], size_u[1]);
    }

    // Print result

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
