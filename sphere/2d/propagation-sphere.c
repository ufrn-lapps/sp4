#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "propagation-sphere.h"
#include "../source.h"

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
    int time_m = 0;
    int time_M = 1000;
    int GRID_SIZE = 100;

    int grid_points[2] = {GRID_SIZE, GRID_SIZE}; // Number of points in each dimension of the grid.
    float points_distance[2] = {1.0, 1.0};       // Distance between each point of the grid in each dimension. in meters.
                                                 // The grid length will be given by grid_points * points_distance

    int x_m = (int)BORDER_SIZE + SPACE_ORDER;
    int x_M = (int)BORDER_SIZE + SPACE_ORDER + GRID_SIZE;
    int y_m = (int)BORDER_SIZE + SPACE_ORDER;
    int y_M = (int)BORDER_SIZE + SPACE_ORDER + GRID_SIZE;

    int size_u[] = {GRID_SIZE + 2 * BORDER_SIZE + 2 * SPACE_ORDER,  // The array u, is composed by the grid points +
                    GRID_SIZE + 2 * BORDER_SIZE + 2 * SPACE_ORDER}; // border on the left and right + space order cells left and right.

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
    float vp_max = 1.5; // In this case, we only have this velocity, but if we have many, we must get the greater.
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

    // Ricker Source
    float fpeak = 1;
    float dt = 1.0;
    float *source = ricker_source(fpeak, time_M - time_m + 1, dt); // dt é o incremento no laço do tempo...
    int source_location[2] = {size_u[0] / 2, size_u[1] / 2};       // Middle_point of each dimension.

    // printf("Fonte: \n");
    // for (int i = 0; i < time_M - time_m + 1; i++)
    // {
    //     printf("%.3f ", source[i]);
    // }
    // printf("\n");

    // printf("Wave propagation...\n");

    float r1 = 0.0784;
    float r0, r2;
    float x_min, x_max;
    float y_min, y_max;
    float raio = 0.0;
    float s1 = source_location[1] * points_distance[1];
    float s2 = source_location[0] * points_distance[0];
    float s3, s4;
    int y_min_coordenate, y_max_coordenate;
    int x_min_coordenate, x_max_coordenate;

    for (int time = time_m, t0 = (time) % (3), t1 = (time + 1) % (3), t2 = (time + 2) % (3);
         time <= time_M;
         time += dt, t0 = (time) % (3), t1 = (time + 1) % (3), t2 = (time + 2) % (3))
    {
        // Calculate dimension extremes
        raio = vp_max * time; // Em metros

        x_min_coordenate = max(y_m + 1, (int)floor((s1 - raio) / points_distance[1]));
        x_max_coordenate = min(y_M + 1, (int)ceil((s1 + raio) / points_distance[1]));

        for (int x = x_min_coordenate; x < x_max_coordenate; x += 1)
        {
            s3 = x - source_location[1] * points_distance[1];
            s4 = sqrt(raio * raio - s3 * s3);

            y_min_coordenate = max(x_m + 1, (int)floor((s2 - s4) / points_distance[0]));
            y_max_coordenate = min(x_M + 1, (int)ceil((s2 + s4) / points_distance[0]));
            // printf("x=(%d,%d), y=(%d, %d)\n", x_min_coordenate, x_max_coordenate, y_min_coordenate, y_max_coordenate);

            for (int y = y_min_coordenate; y < y_max_coordenate; y += 1)
            {
                r0 = vp[x][y] * vp[x][y];
                r2 = r0 * r1;
                u[t1][x][y] = -4.0F * r2 * u[t0][x][y] +
                              r2 * (u[t0][x - 1][y] +
                                    u[t0][x][y - 1] +
                                    u[t0][x][y + 1] +
                                    u[t0][x + 1][y] -
                                    u[t2][x][y]) +
                              2.0F * u[t0][x][y];
            }
        }
        // Inject source
        u[t1][source_location[0]][source_location[1]] -= source[time - 1] * vp[source_location[0]][source_location[1]] * vp[source_location[0]][source_location[1]] * dt * dt;

        // if (time == 53)
        // {
        //     printf("x=(%d, %d), y=(%d, %d)\n", x_min_coordenate, x_max_coordenate, y_min_coordenate, y_max_coordenate);
        //     print_array_2d(u[t1], size_u[0], size_u[1]);
        //     exit(0);
        // }
        // if (y_min == (y_m - 1) && y_max == (y_M - 1) &&
        //     x_min == (x_m - 1) && x_max == (x_M - 1) && alcancou == 0)
        // {
        //     printf("Na iteração %d foi alncaçando os limites da grid.\n", time);
        //     alcancou = 1;
        // }

        // if (time == 1 || time == 10 || time == 100 || time == 1000 || time == 10000)
        // {
        //     printf("Numero de pontos calculados no tempo=%d: %d\n", time, count);
        // }
    }

    // Print result
    print_array_2d(u[1], size_u[0], size_u[1]);

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
