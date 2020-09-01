#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "propagation-reference-3d.h"
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
int SPACE_ORDER = 2;

int main(int argc, char const *argv[])
{
    int BORDER_SIZE = 40;
    int time_m = 0;
    int time_M = 1000;
    int GRID_SIZE = 100;

    int x_m = 0;
    int x_M = 2 * BORDER_SIZE + GRID_SIZE;
    int y_m = 0;
    int y_M = 2 * BORDER_SIZE + GRID_SIZE;
    int z_m = 0;
    int z_M = 2 * BORDER_SIZE + GRID_SIZE;

    int u_size[] = {GRID_SIZE + 2 * BORDER_SIZE + 2 * SPACE_ORDER + 1,  // The array u, is composed by the grid points +
                    GRID_SIZE + 2 * BORDER_SIZE + 2 * SPACE_ORDER + 1,  // border on the left and right + space order cells left and right.
                    GRID_SIZE + 2 * BORDER_SIZE + 2 * SPACE_ORDER + 1}; // O +1 é pq o loop do espaço é menor igual. 100 + 2*40 + 2*2 + 1 = 185
    float ***vp;
    float ***damp;
    float ****u;

    // Allocate data
    // printf("Allocating u array with dimensions %d x %d x %d...\n", TIME_ORDER + 1, u_size[0], u_size[1]);
    u = (float ****)malloc(sizeof(float *) * (TIME_ORDER + 1));
    for (int i = 0; i < TIME_ORDER + 1; i++)
    {
        u[i] = (float ***)malloc(sizeof(float *) * u_size[0]);
        for (int j = 0; j < u_size[0]; j++)
        {
            u[i][j] = (float **)malloc(sizeof(float *) * u_size[1]);
            for (int k = 0; k < u_size[1]; k++)
            {
                u[i][j][k] = (float *)malloc(sizeof(float) * u_size[2]);
            }
        }
    }

    // printf("Allocating vp array with dimensions %d x %d...\n", u_size[0], u_size[1]);
    float vp_max = 1.5; // In this case, we only have this velocity, but if we have many, we must get the greater.

    vp = (float ***)malloc(sizeof(float **) * u_size[0]);
    for (int i; i < u_size[0]; i++)
    {
        vp[i] = (float **)malloc(sizeof(float **) * u_size[1]);
        for (int j = 0; j < u_size[1]; j++)
        {
            vp[i][j] = (float *)malloc(sizeof(float) * u_size[2]);
        }
    }

    // Initializing values
    // printf("Initializing data...\n");
    for (int i = 0; i < u_size[0]; i++)
    {
        for (int j = 0; j < u_size[1]; j++)
        {
            for (int k = 0; k < u_size[2]; k++)
            {
                u[0][i][j][k] = 0.0;
                u[1][i][j][k] = 0.0;
                u[2][i][j][k] = 0.0;
                vp[i][j][k] = vp_max;
            }
        }
    }

    // Damp
    damp = initialize_damp(BORDER_SIZE,
                           u_size[0] - 2 * SPACE_ORDER,
                           u_size[1] - 2 * SPACE_ORDER,
                           u_size[2] - 2 * SPACE_ORDER);
    // print_array_3d(damp, u_size[0], u_size[1], u_size[2]);

    // Ricker Source
    float fpeak = 1;
    float dt = 0.152;
    float *source = ricker_source(fpeak, time_M - time_m + 1, dt);
    int source_location[3] = {(int)floor(u_size[0] / 2),
                              (int)floor(u_size[1] / 2),
                              (int)floor(u_size[2] / 2)};

    float r1 = dt * dt;
    float r0;
    // Adicionado para a injeção da fonte.
    float o_x = -40.0;
    float o_y = -40.0;
    float o_z = -40.0;

    int time, t0, t1, t2;

    for (time = time_m, t0 = (time) % (3), t1 = (time + 1) % (3), t2 = (time + 2) % (3);
         time <= time_M;
         time += 1, t0 = (time) % (3), t1 = (time + 1) % (3), t2 = (time + 2) % (3))
    {
        for (int x = x_m; x <= x_M; x += 1)
        {
            for (int y = y_m; y <= y_M; y += 1)
            {
                for (int z = z_m; z <= z_M; z += 1)
                {
                    r0 = vp[x + 2][y + 2][z + 2] * vp[x + 2][y + 2][z + 2];
                    u[t1][x + 2][y + 2][z + 2] = 2.0F * (r0 * (-3.0F * r1 * u[t0][x + 2][y + 2][z + 2] + 5.0e-1F * (r1 * u[t0][x + 1][y + 2][z + 2] + r1 * u[t0][x + 2][y + 1][z + 2] + r1 * u[t0][x + 2][y + 2][z + 1] + r1 * u[t0][x + 2][y + 2][z + 3] + r1 * u[t0][x + 2][y + 3][z + 2] + r1 * u[t0][x + 3][y + 2][z + 2] + dt * damp[x + 1][y + 1][z + 1] * u[t0][x + 2][y + 2][z + 2])) + 1.0F * u[t0][x + 2][y + 2][z + 2] - 5.0e-1F * u[t2][x + 2][y + 2][z + 2]) / (r0 * dt * damp[x + 1][y + 1][z + 1] + 1);
                }
            }
        }

        int ii_src_0 = (int)(floor(-1.0 * o_x + 1.0 * source_location[0]));
        int ii_src_1 = (int)(floor(-1.0 * o_y + 1.0 * source_location[1]));
        int ii_src_2 = (int)(floor(-1.0 * o_z + 1.0 * source_location[2]));
        int ii_src_3 = (int)(floor(-1.0 * o_z + 1.0 * source_location[2])) + 1;
        int ii_src_4 = (int)(floor(-1.0 * o_y + 1.0 * source_location[1])) + 1;
        int ii_src_5 = (int)(floor(-1.0 * o_x + 1.0 * source_location[0])) + 1;

        float px = (float)(-o_x - 1.0F * (int)(floor(-1.0F * o_x + 1.0F * source_location[0])) + source_location[0]);
        float py = (float)(-o_y - 1.0F * (int)(floor(-1.0F * o_y + 1.0F * source_location[1])) + source_location[1]);
        float pz = (float)(-o_z - 1.0F * (int)(floor(-1.0F * o_z + 1.0F * source_location[2])) + source_location[2]);

        if (ii_src_0 >= x_m - 1 && ii_src_1 >= y_m - 1 && ii_src_2 >= z_m - 1 && ii_src_0 <= x_M + 1 && ii_src_1 <= y_M + 1 && ii_src_2 <= z_M + 1)
        {
            float r2 = 2.31039985141754e-2F * (vp[ii_src_0 + 2][ii_src_1 + 2][ii_src_2 + 2] * vp[ii_src_0 + 2][ii_src_1 + 2][ii_src_2 + 2]) * (-1.0F * px * py * pz + 1.0F * px * py + 1.0F * px * pz - 1.0F * px + 1.0F * py * pz - 1.0F * py - 1.0F * pz + 1) * source[time];
            u[t1][ii_src_0 + 2][ii_src_1 + 2][ii_src_2 + 2] += r2;
            // printf("Injecting -> u[%d][%d][%d][%d]=%f\n", time, ii_src_0 + 2, ii_src_1 + 2, ii_src_2 + 2, r2);
        }
        if (ii_src_0 >= x_m - 1 && ii_src_1 >= y_m - 1 && ii_src_3 >= z_m - 1 && ii_src_0 <= x_M + 1 && ii_src_1 <= y_M + 1 && ii_src_3 <= z_M + 1)
        {
            float r3 = 2.31039985141754e-2F * (vp[ii_src_0 + 2][ii_src_1 + 2][ii_src_3 + 2] * vp[ii_src_0 + 2][ii_src_1 + 2][ii_src_3 + 2]) * (1.0F * px * py * pz - 1.0F * px * pz - 1.0F * py * pz + 1.0F * pz) * source[time];
            u[t1][ii_src_0 + 2][ii_src_1 + 2][ii_src_3 + 2] += r3;
        }
        if (ii_src_0 >= x_m - 1 && ii_src_2 >= z_m - 1 && ii_src_4 >= y_m - 1 && ii_src_0 <= x_M + 1 && ii_src_2 <= z_M + 1 && ii_src_4 <= y_M + 1)
        {
            float r4 = 2.31039985141754e-2F * (vp[ii_src_0 + 2][ii_src_4 + 2][ii_src_2 + 2] * vp[ii_src_0 + 2][ii_src_4 + 2][ii_src_2 + 2]) * (1.0F * px * py * pz - 1.0F * px * py - 1.0F * py * pz + 1.0F * py) * source[time];
            u[t1][ii_src_0 + 2][ii_src_4 + 2][ii_src_2 + 2] += r4;
        }
        if (ii_src_0 >= x_m - 1 && ii_src_3 >= z_m - 1 && ii_src_4 >= y_m - 1 && ii_src_0 <= x_M + 1 && ii_src_3 <= z_M + 1 && ii_src_4 <= y_M + 1)
        {
            float r5 = 2.31039985141754e-2F * (vp[ii_src_0 + 2][ii_src_4 + 2][ii_src_3 + 2] * vp[ii_src_0 + 2][ii_src_4 + 2][ii_src_3 + 2]) * (-1.0F * px * py * pz + 1.0F * py * pz) * source[time];
            u[t1][ii_src_0 + 2][ii_src_4 + 2][ii_src_3 + 2] += r5;
        }
        if (ii_src_1 >= y_m - 1 && ii_src_2 >= z_m - 1 && ii_src_5 >= x_m - 1 && ii_src_1 <= y_M + 1 && ii_src_2 <= z_M + 1 && ii_src_5 <= x_M + 1)
        {
            float r6 = 2.31039985141754e-2F * (vp[ii_src_5 + 2][ii_src_1 + 2][ii_src_2 + 2] * vp[ii_src_5 + 2][ii_src_1 + 2][ii_src_2 + 2]) * (1.0F * px * py * pz - 1.0F * px * py - 1.0F * px * pz + 1.0F * px) * source[time];
            u[t1][ii_src_5 + 2][ii_src_1 + 2][ii_src_2 + 2] += r6;
        }
        if (ii_src_1 >= y_m - 1 && ii_src_3 >= z_m - 1 && ii_src_5 >= x_m - 1 && ii_src_1 <= y_M + 1 && ii_src_3 <= z_M + 1 && ii_src_5 <= x_M + 1)
        {
            float r7 = 2.31039985141754e-2F * (vp[ii_src_5 + 2][ii_src_1 + 2][ii_src_3 + 2] * vp[ii_src_5 + 2][ii_src_1 + 2][ii_src_3 + 2]) * (-1.0F * px * py * pz + 1.0F * px * pz) * source[time];
            u[t1][ii_src_5 + 2][ii_src_1 + 2][ii_src_3 + 2] += r7;
        }
        if (ii_src_2 >= z_m - 1 && ii_src_4 >= y_m - 1 && ii_src_5 >= x_m - 1 && ii_src_2 <= z_M + 1 && ii_src_4 <= y_M + 1 && ii_src_5 <= x_M + 1)
        {
            float r8 = 2.31039985141754e-2F * (vp[ii_src_5 + 2][ii_src_4 + 2][ii_src_2 + 2] * vp[ii_src_5 + 2][ii_src_4 + 2][ii_src_2 + 2]) * (-1.0F * px * py * pz + 1.0F * px * py) * source[time];
            u[t1][ii_src_5 + 2][ii_src_4 + 2][ii_src_2 + 2] += r8;
        }
        if (ii_src_3 >= z_m - 1 && ii_src_4 >= y_m - 1 && ii_src_5 >= x_m - 1 && ii_src_3 <= z_M + 1 && ii_src_4 <= y_M + 1 && ii_src_5 <= x_M + 1)
        {
            float r9 = 2.31039985141754e-2F * px * py * pz * (vp[ii_src_5 + 2][ii_src_4 + 2][ii_src_3 + 2] * vp[ii_src_5 + 2][ii_src_4 + 2][ii_src_3 + 2]) * source[time];
            u[t1][ii_src_5 + 2][ii_src_4 + 2][ii_src_3 + 2] += r9;
        }
    }

    print_array_3d(u[t2], u_size[0], u_size[1], u_size[2]);

    // Free resources
    free(source);
    free(vp);
    free(u[0]);
    free(u[1]);
    free(u[2]);

    return 0;
}

void print_array_3d(float ***u, int x_size, int y_size, int z_size)
{

    for (int i = 0; i < x_size; i++)
    {
        for (int j = 0; j < y_size; j++)
        {
            for (int k = 0; k < z_size; k++)
            {
                printf("%g ", u[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
    printf("\n\n");
}

float ***initialize_damp(int border_size, int size_x, int size_y, int size_z)
{
    float ***damp;

    damp = (float ***)malloc(sizeof(float **) * size_x);
    for (int i = 0; i <= size_x; i++)
    {
        damp[i] = (float **)malloc(sizeof(float *) * size_y);
        for (int j = 0; j <= size_y; j++)
        {
            damp[i][j] = (float *)malloc(sizeof(float) * size_z);
        }
    }

    int x_m = 0;
    int x_M = 180;
    int y_m = 0;
    int y_M = 180;
    int z_m = 0;
    int z_M = 180;
    float h_x = 1.0;
    float h_y = 1.0;
    float h_z = 1.0;
    int abc_z_r_rtkn = 40;
    int abc_z_l_ltkn = 40;
    int abc_y_r_rtkn = 40;
    int abc_y_l_ltkn = 40;
    int abc_x_r_rtkn = 40;
    int abc_x_l_ltkn = 40;

    for (int abc_x_l = x_m; abc_x_l <= abc_x_l_ltkn + x_m - 1; abc_x_l += 1)
    {
        for (int y = y_m; y <= y_M; y += 1)
        {
            for (int z = z_m; z <= z_M; z += 1)
            {
                damp[abc_x_l + 1][y + 1][z + 1] += (-4.12276274369678e-2F * sin(6.28318530717959F * fabs(2.5e-2F * x_m - 2.5e-2F * abc_x_l + 1.025F)) + 2.5904082296183e-1F * fabs(2.5e-2F * x_m - 2.5e-2F * abc_x_l + 1.025F)) / h_x;
            }
        }
    }

    for (int abc_x_r = -abc_x_r_rtkn + x_M + 1; abc_x_r <= x_M; abc_x_r += 1)
    {
        for (int y = y_m; y <= y_M; y += 1)
        {
            for (int z = z_m; z <= z_M; z += 1)
            {
                damp[abc_x_r + 1][y + 1][z + 1] += (-4.12276274369678e-2F * sin(6.28318530717959F * fabs(-2.5e-2F * x_M + 2.5e-2F * abc_x_r + 1.025F)) + 2.5904082296183e-1F * fabs(-2.5e-2F * x_M + 2.5e-2F * abc_x_r + 1.025F)) / h_x;
            }
        }
    }
    for (int x = x_m; x <= x_M; x += 1)
    {
        for (int abc_y_l = y_m; abc_y_l <= abc_y_l_ltkn + y_m - 1; abc_y_l += 1)
        {
            for (int z = z_m; z <= z_M; z += 1)
            {
                damp[x + 1][abc_y_l + 1][z + 1] += (-4.12276274369678e-2F * sin(6.28318530717959F * fabs(2.5e-2F * y_m - 2.5e-2F * abc_y_l + 1.025F)) + 2.5904082296183e-1F * fabs(2.5e-2F * y_m - 2.5e-2F * abc_y_l + 1.025F)) / h_y;
            }
        }
        for (int abc_y_r = -abc_y_r_rtkn + y_M + 1; abc_y_r <= y_M; abc_y_r += 1)
        {
            for (int z = z_m; z <= z_M; z += 1)
            {
                damp[x + 1][abc_y_r + 1][z + 1] += (-4.12276274369678e-2F * sin(6.28318530717959F * fabs(-2.5e-2F * y_M + 2.5e-2F * abc_y_r + 1.025F)) + 2.5904082296183e-1F * fabs(-2.5e-2F * y_M + 2.5e-2F * abc_y_r + 1.025F)) / h_y;
            }
        }
        for (int y = y_m; y <= y_M; y += 1)
        {
            for (int abc_z_l = z_m; abc_z_l <= abc_z_l_ltkn + z_m - 1; abc_z_l += 1)
            {
                damp[x + 1][y + 1][abc_z_l + 1] += (-4.12276274369678e-2F * sin(6.28318530717959F * fabs(2.5e-2F * z_m - 2.5e-2F * abc_z_l + 1.025F)) + 2.5904082296183e-1F * fabs(2.5e-2F * z_m - 2.5e-2F * abc_z_l + 1.025F)) / h_z;
            }
            for (int abc_z_r = -abc_z_r_rtkn + z_M + 1; abc_z_r <= z_M; abc_z_r += 1)
            {
                damp[x + 1][y + 1][abc_z_r + 1] += (-4.12276274369678e-2F * sin(6.28318530717959F * fabs(-2.5e-2F * z_M + 2.5e-2F * abc_z_r + 1.025F)) + 2.5904082296183e-1F * fabs(-2.5e-2F * z_M + 2.5e-2F * abc_z_r + 1.025F)) / h_z;
            }
        }
    }

    return damp;
}
