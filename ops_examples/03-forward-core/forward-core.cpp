#define _POSIX_C_SOURCE 200809L
#include "math.h"
// #include "pmmintrin.h"
#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"
// #include "xmmintrin.h"

struct dataobj {
    void *restrict data;
    int *size;
    int *npsize;
    int *dsize;
    int *hsize;
    int *hofs;
    int *oofs;
};

struct profiler {
    double section0;
    double section1;
    double section2;
};

void print_to_file(char *filename, float *mydata, int size_x, int size_y);

int Forward(struct dataobj *restrict damp_vec, const float dt, struct dataobj *restrict m_vec, const float o_x,
            const float o_y, struct dataobj *restrict rec_vec, struct dataobj *restrict rec_coords_vec,
            struct dataobj *restrict src_vec, struct dataobj *restrict src_coords_vec, struct dataobj *restrict u_vec,
            const int x_M, const int x_m, const int y_M, const int y_m, const int p_rec_M, const int p_rec_m,
            const int p_src_M, const int p_src_m, const int time_M, const int time_m, struct profiler *timers) {
    float(*restrict damp)[damp_vec->size[1]] __attribute__((aligned(64))) = (float(*)[damp_vec->size[1]])damp_vec->data;
    float(*restrict m)[m_vec->size[1]] __attribute__((aligned(64))) = (float(*)[m_vec->size[1]])m_vec->data;
    float(*restrict rec)[rec_vec->size[1]] __attribute__((aligned(64))) = (float(*)[rec_vec->size[1]])rec_vec->data;
    float(*restrict rec_coords)[rec_coords_vec->size[1]] __attribute__((aligned(64))) =
        (float(*)[rec_coords_vec->size[1]])rec_coords_vec->data;
    float(*restrict src)[src_vec->size[1]] __attribute__((aligned(64))) = (float(*)[src_vec->size[1]])src_vec->data;
    float(*restrict src_coords)[src_coords_vec->size[1]] __attribute__((aligned(64))) =
        (float(*)[src_coords_vec->size[1]])src_coords_vec->data;
    float(*restrict u)[u_vec->size[1]][u_vec->size[2]] __attribute__((aligned(64))) =
        (float(*)[u_vec->size[1]][u_vec->size[2]])u_vec->data;

    // Print input to file.
    print_to_file("input-rec0/damp.txt", damp, damp_vec->size[1], damp_vec->size[0]);
    print_to_file("input-rec0/m.txt", m, m_vec->size[1], m_vec->size[0]);
    print_to_file("input-rec0/rec.txt", rec, rec_vec->size[1], rec_vec->size[0]);
    print_to_file("input-rec0/rec_coords.txt", rec_coords, rec_coords_vec->size[1], rec_coords_vec->size[0]);
    print_to_file("input-rec0/src.txt", src, src_vec->size[1], src_vec->size[0]);
    print_to_file("input-rec0/src_coords.txt", src_coords, src_coords_vec->size[1], src_coords_vec->size[0]);
    // print_to_file("input/u_vec.txt", u, u_vec);

    /* Flush denormal numbers to zero in hardware */
    // _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    // _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    // printf("time_m=%d time_M=%d\n", time_m, time_M);
    // printf("x_m=%d  x_M=%d y_m=%d y_M=%d\n", x_m, x_M, y_m, y_M);
    // printf("p_rec_m=%d p_rec_M=%d\n", p_rec_m, p_rec_M);
    for (int time = time_m; time <= time_M; time += 1) {
        // printf("------------------------------------ T = %d -----------------------------------\n", time);
        struct timeval start_section0, end_section0;
        gettimeofday(&start_section0, NULL);
        /* Begin section0 */
        for (int x = x_m; x <= x_M; x += 1) {
#pragma omp simd aligned(damp, m, u : 32)
            for (int y = y_m; y <= y_M; y += 1) {
                u[time + 1][x + 4][y + 4] =
                    2.0F *
                    (2.5e-1F * dt * damp[x + 1][y + 1] * u[time - 1][x + 4][y + 4] +
                     (dt * dt) * (-4.16666657e-4F * (u[time][x + 2][y + 4] + u[time][x + 4][y + 2] +
                                                     u[time][x + 4][y + 6] + u[time][x + 6][y + 4]) +
                                  6.66666652e-3F * (u[time][x + 3][y + 4] + u[time][x + 4][y + 3] +
                                                    u[time][x + 4][y + 5] + u[time][x + 5][y + 4]) -
                                  2.49999994e-2F * u[time][x + 4][y + 4]) +
                     1.0F * m[x + 2][y + 2] * u[time][x + 4][y + 4] -
                     5.0e-1F * m[x + 2][y + 2] * u[time - 1][x + 4][y + 4]) /
                    (5.0e-1F * dt * damp[x + 1][y + 1] + 1.0F * m[x + 2][y + 2]);

                // printf("u[%d][%d]=%f\n", x + 4, y + 4, u[time + 1][x + 4][y + 4]);

                // if ((x == 80 && y == 90) || (x == 80 && y == 78))
                // if (x == 2)
                // {
                //     printf("damp=%f\tu[%d][%d]=%e\n", damp[x + 1][y + 1], x + 4, y + 4, u[time + 1][x + 4][y + 4]);
                // }
            }
        }
        /* End section0 */
        gettimeofday(&end_section0, NULL);
        timers->section0 += (double)(end_section0.tv_sec - start_section0.tv_sec) +
                            (double)(end_section0.tv_usec - start_section0.tv_usec) / 1000000;
        struct timeval start_section1, end_section1;
        gettimeofday(&start_section1, NULL);
        /* Begin section1 */
        // printf("p_src_m=%d p_src_M=%d\n", p_src_m, p_src_M);
        for (int p_src = p_src_m; p_src <= p_src_M; p_src += 1) {
            float r0 = (int)(floor(-1.0e-1F * o_x + 1.0e-1F * src_coords[p_src][0]));
            int ii_src_0 = r0;
            float r1 = (int)(floor(-1.0e-1F * o_y + 1.0e-1F * src_coords[p_src][1]));
            int ii_src_1 = r1;
            int ii_src_2 = r1 + 1;
            int ii_src_3 = r0 + 1;

            float px = (float)(-1.0e+1F * r0 - o_x + src_coords[p_src][0]);
            float py = (float)(-1.0e+1F * r1 - o_y + src_coords[p_src][1]);
            // printf("r0=%f r1=%f - o_x=%f o_y=%f - src_coords=%f %f\n", r0, r1, o_x, o_y, src_coords[p_src][0],
            //        src_coords[p_src][1]);
            // printf("p[0]=%f p[1]-%f\n", py, px);
            // printf("x_m=%d, x_M=%d, y_m=%d, y_M=%d\n", x_m, x_M, y_m, y_M);

            // printf("p_src=%d\tsrc_coords=%f\tsrc=%f\n", p_src, src_coords[p_src][0], src[time][p_src]);
            if (ii_src_0 >= x_m - 1 && ii_src_1 >= y_m - 1 && ii_src_0 <= x_M + 1 && ii_src_1 <= y_M + 1) {
                int r2 = ii_src_0 + 4;
                int r3 = ii_src_1 + 4;
                int r4 = ii_src_0 + 2;
                int r5 = ii_src_1 + 2;
                float r6 =
                    (dt * dt) * (1.0e-2F * px * py - 1.0e-1F * px - 1.0e-1F * py + 1) * src[time][p_src] / m[r4][r5];
                u[time + 1][r2][r3] += r6;

                // if (time == 1) {
                // printf("dt=%f px=%f py=%f src_value=%f m=%f\n", dt, px, py, src[time][p_src], m[r4][r5]);
                // printf("r2=%d r3=%d\n", r2, r3);
                // printf("u[%d][%d] injected=%f\n", r2, r3, r6);
                // }
            }
            if (ii_src_0 >= x_m - 1 && ii_src_2 >= y_m - 1 && ii_src_0 <= x_M + 1 && ii_src_2 <= y_M + 1) {
                int r7 = ii_src_0 + 4;
                int r8 = ii_src_2 + 4;
                int r9 = ii_src_0 + 2;
                int r10 = ii_src_2 + 2;
                float r11 = (dt * dt) * (-1.0e-2F * px * py + 1.0e-1F * py) * src[time][p_src] / m[r9][r10];
                u[time + 1][r7][r8] += r11;
                // printf("u[%d][%d] injected=%f\n", r7, r8, r11);
            }
            if (ii_src_1 >= y_m - 1 && ii_src_3 >= x_m - 1 && ii_src_1 <= y_M + 1 && ii_src_3 <= x_M + 1) {
                int r12 = ii_src_3 + 4;
                int r13 = ii_src_1 + 4;
                int r14 = ii_src_3 + 2;
                int r15 = ii_src_1 + 2;
                float r16 = (dt * dt) * (-1.0e-2F * px * py + 1.0e-1F * px) * src[time][p_src] / m[r14][r15];
                // printf("r12=%d r13=%d\n", r12, r13);
                u[time + 1][r12][r13] += r16;
                // printf("u[%d][%d] injected=%f\n", r12, r13, r16);
            }
            if (ii_src_2 >= y_m - 1 && ii_src_3 >= x_m - 1 && ii_src_2 <= y_M + 1 && ii_src_3 <= x_M + 1) {
                int r17 = ii_src_3 + 4;
                int r18 = ii_src_2 + 4;
                int r19 = ii_src_3 + 2;
                int r20 = ii_src_2 + 2;
                float r21 = 1.0e-2F * px * py * (dt * dt) * src[time][p_src] / m[r19][r20];
                // printf("r17=%d r18=%d\n", r17, r18);
                u[time + 1][r17][r18] += r21;
                // printf("u[%d][%d] injected=%f\n", r17, r18, r21);
            }
        }

        // if (time == 2) print_to_file("output-core/u-partial-2.txt", u[time + 1], u_vec->size[2], u_vec->size[1]);
        // if (time == 10) print_to_file("output-core/u-partial-10.txt", u[time + 1], u_vec->size[2], u_vec->size[1]);
        // if (time == 50) print_to_file("output-core/u-partial-50.txt", u[time + 1], u_vec->size[2], u_vec->size[1]);
        // if (time == 100) print_to_file("output-core/u-partial-100.txt", u[time + 1], u_vec->size[2], u_vec->size[1]);

        /* End section1 */
        gettimeofday(&end_section1, NULL);
        timers->section1 += (double)(end_section1.tv_sec - start_section1.tv_sec) +
                            (double)(end_section1.tv_usec - start_section1.tv_usec) / 1000000;
        struct timeval start_section2, end_section2;
        gettimeofday(&start_section2, NULL);
        /* Begin section2 */
        for (int p_rec = p_rec_m; p_rec <= p_rec_M; p_rec += 1) {
            float r22 = (int)(floor(-1.0e-1F * o_x + 1.0e-1F * rec_coords[p_rec][0]));
            int ii_rec_0 = r22;
            float r23 = (int)(floor(-1.0e-1F * o_y + 1.0e-1F * rec_coords[p_rec][1]));
            int ii_rec_1 = r23;
            int ii_rec_2 = r23 + 1;
            int ii_rec_3 = r22 + 1;

            float px = (float)(-1.0e+1F * r22 - o_x + rec_coords[p_rec][0]);
            float py = (float)(-1.0e+1F * r23 - o_y + rec_coords[p_rec][1]);

            // printf("ii_rec=[%d %d %d %d] p=[%f %f]\n", ii_rec_0, ii_rec_1, ii_rec_2, ii_rec_3, px, py);

            float sum = 0.0F;
            if (ii_rec_0 >= x_m - 1 && ii_rec_1 >= y_m - 1 && ii_rec_0 <= x_M + 1 && ii_rec_1 <= y_M + 1) {
                int r24 = ii_rec_0 + 4;
                int r25 = ii_rec_1 + 4;
                sum += (1.0e-2F * px * py - 1.0e-1F * px - 1.0e-1F * py + 1) * u[time][r24][r25];
                // printf("u = %e p[0] = %f  p[1] = %f\n", u[time][r24][r25], px, py);
            }
            if (ii_rec_0 >= x_m - 1 && ii_rec_2 >= y_m - 1 && ii_rec_0 <= x_M + 1 && ii_rec_2 <= y_M + 1) {
                int r26 = ii_rec_0 + 4;
                int r27 = ii_rec_2 + 4;
                sum += (-1.0e-2F * px * py + 1.0e-1F * py) * u[time][r26][r27];
            }
            if (ii_rec_1 >= y_m - 1 && ii_rec_3 >= x_m - 1 && ii_rec_1 <= y_M + 1 && ii_rec_3 <= x_M + 1) {
                int r28 = ii_rec_3 + 4;
                int r29 = ii_rec_1 + 4;
                sum += (-1.0e-2F * px * py + 1.0e-1F * px) * u[time][r28][r29];
            }
            if (ii_rec_2 >= y_m - 1 && ii_rec_3 >= x_m - 1 && ii_rec_2 <= y_M + 1 && ii_rec_3 <= x_M + 1) {
                int r30 = ii_rec_3 + 4;
                int r31 = ii_rec_2 + 4;
                sum += 1.0e-2F * px * py * u[time][r30][r31];
            }
            rec[time][p_rec] = sum;
            // printf("rec position %d x %d value = %e \t", time, p_rec, sum);
        }
        // printf("\n");
        /* End section2 */
        gettimeofday(&end_section2, NULL);
        timers->section2 += (double)(end_section2.tv_sec - start_section2.tv_sec) +
                            (double)(end_section2.tv_usec - start_section2.tv_usec) / 1000000;
    }

    // printf("rec core size: %d x %d x %d\n", rec_vec->size[2], rec_vec->size[1], rec_vec->size[0]);
    print_to_file("output-core/rec.txt", rec, rec_vec->size[1], rec_vec->size[0]);
    print_to_file("output-core/u.txt", u[time_M + 1], u_vec->size[2], u_vec->size[1]);

    return 0;
}

/**
 * @brief  Write data structure into a file.
 * @note   If file exist it will be overwritten.
 * @param  *filename: Name of the file to be saved.
 * @param  *mydata: Data to be written in the file. Float array expected.
 * @param  metadata: struct with size of the data.
 * @retval None
 */
void print_to_file(char *filename, float *mydata, int size_x, int size_y) {
    FILE *fptr;
    fptr = fopen(filename, "w");

    if (fptr == NULL) {
        printf("Error!");
        exit(1);
    }

    // Write meta data into the file.
    printf("Filename=%s\tSize=%dx%d\n", filename, size_x, size_y);
    fprintf(fptr, "%d %d\n", size_x, size_y);

    // Write data.
    for (int j = 0; j < size_y; j++) {
        for (int i = 0; i < size_x; i++) {
            fprintf(fptr, "%e ", mydata[j * size_x + i]);
        }
        fprintf(fptr, "\n");
    }

    fclose(fptr);
    return;
}