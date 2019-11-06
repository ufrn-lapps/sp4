#define _POSIX_C_SOURCE 200809L
#include "math.h"
#include "pmmintrin.h"
#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"
#include "xmmintrin.h"

#define OPS_2D
#include <ops_seq.h>

#include "forward-ops.h"
#include "utils.h"

// Define global variables here
int SPACE_ORDER = 4;
int GRID_SIZE = 101;
int BORDER_SIZE = 40;
int T_INTERVALS = 714;

int main(int argc, char *argv[]) {
    // Input variables
    float *damp, *m, *src, *rec, *src_coords, *rec_coords, **u;
    float dt;
    int size_u[] = {GRID_SIZE + 2 * BORDER_SIZE + 2 * SPACE_ORDER, GRID_SIZE + 2 * BORDER_SIZE + 2 * SPACE_ORDER};
    int size_m[] = {GRID_SIZE + 2 * BORDER_SIZE + SPACE_ORDER, GRID_SIZE + 2 * BORDER_SIZE + SPACE_ORDER};
    int size_damp[] = {GRID_SIZE + 2 * BORDER_SIZE + 2, GRID_SIZE + 2 * BORDER_SIZE + 2};
    int size_rec[] = {GRID_SIZE, GRID_SIZE};
    int ii_src[4], ii_rec[4];  // Source position, and source interpolated positions
    float p[2];

    /*****************************************************************************************************************/
    /**************                                   Read input variables                               *************/
    /*****************************************************************************************************************/

    damp = read_2d("input-rec0/damp.txt");
    m = read_2d("input-rec0/m.txt");
    src = read_2d("input-rec0/src.txt");
    rec = read_2d("input-rec0/rec.txt");
    src_coords = read_2d("input-rec0/src_coords.txt");
    rec_coords = read_2d("input-rec0/rec_coords.txt");

    float o_x = -400.00, o_y = -400.00;
    int x_limits[] = {0, 180};
    int y_limits[] = {0, 180};

    int p_src_m = 0;
    int p_src_M = 0;
    int p_rec_m = 0;
    int p_rec_M = 100;

    u = (float **)malloc(sizeof(float *) * 3);
    u[0] = (float *)malloc(sizeof(float) * size_u[0] * size_u[1]);
    u[1] = (float *)malloc(sizeof(float) * size_u[0] * size_u[1]);
    u[2] = (float *)malloc(sizeof(float) * size_u[0] * size_u[1]);

    dt = 1.4;

    /*****************************************************************************************************************/
    /**************                                 Initialize OPS variables                             *************/
    /*****************************************************************************************************************/

    ops_init(0, NULL, 3);

    // Declare ops_block
    ops_block grid = ops_decl_block(2, "grid");

    // Declare ops_dat objects
    int base[] = {0, 0};
    int d_p[] = {0, 0};
    int d_m[] = {0, 0};

    ops_dat dat_ut[3];
    ops_dat dat_m;
    ops_dat dat_damp;
    ops_reduction red_sum;

    dat_ut[0] = ops_decl_dat(grid, 1, size_u, base, d_m, d_p, u[0], "float", "ut0");
    dat_ut[1] = ops_decl_dat(grid, 1, size_u, base, d_m, d_p, u[1], "float", "ut1");
    dat_ut[2] = ops_decl_dat(grid, 1, size_u, base, d_m, d_p, u[2], "float", "ut2");
    dat_m = ops_decl_dat(grid, 1, size_m, base, d_m, d_p, m, "float", "m");
    dat_damp = ops_decl_dat(grid, 1, size_damp, base, d_m, d_p, damp, "float", "damp");
    red_sum = ops_decl_reduction_handle(sizeof(float), "float", "sum");

    // Declare stencils
    int stencil_1_point[] = {0, 0};
    int stencil_1F_point[] = {1, 1};
    int stencil_m2_m2[] = {-2, -2};
    int stencil_m3_m3[] = {-3, -3};
    int stencil_9_points[] = {0, 1, 0, 0, -1, 0, -2, 0, 2, 0, 0, -2, 0, -1, 1, 0, 0, 2};

    ops_stencil S2D_00 = ops_decl_stencil(2, 1, stencil_1_point, "0,0");
    ops_stencil S2D_m2_m2 = ops_decl_stencil(2, 1, stencil_m2_m2, "-2,-2");
    ops_stencil S2D_m3_m3 = ops_decl_stencil(2, 1, stencil_m3_m3, "-3,-3");
    ops_stencil S2D_11 = ops_decl_stencil(2, 1, stencil_1F_point, "1,1");
    ops_stencil S2D_SO4 = ops_decl_stencil(2, 9, stencil_9_points, "so4");

    // Define kernel ranges
    int propagation_range[] = {SPACE_ORDER, SPACE_ORDER + GRID_SIZE + 2 * BORDER_SIZE, SPACE_ORDER,
                               SPACE_ORDER + GRID_SIZE + 2 * BORDER_SIZE};

    ops_partition("");

    // printf("x_m=%d x_M=%d y_m=%d y_M=%d\n", SPACE_ORDER, SPACE_ORDER + GRID_SIZE + 2 * BORDER_SIZE, SPACE_ORDER,
    //        SPACE_ORDER + GRID_SIZE + 2 * BORDER_SIZE);
    // printf("time_m=%d time_M=%d\n", 1, T_INTERVALS);
    // ops_printf("p_rec_m = %d, p_rec_M=%d\n", p_rec_m, p_rec_M);
    for (int t = 1; t <= T_INTERVALS; t += 1) {
        // printf("------------------------------------ T = %d -----------------------------------\n", t);
        ops_par_loop(propagation_kernel, "propagation_kernel", grid, 2, propagation_range,
                     ops_arg_dat(dat_m, 1, S2D_m2_m2, "float", OPS_READ),
                     ops_arg_dat(dat_damp, 1, S2D_m3_m3, "float", OPS_READ),
                     ops_arg_dat(dat_ut[t % 3], 1, S2D_SO4, "float", OPS_READ),        // t
                     ops_arg_dat(dat_ut[(t - 1) % 3], 1, S2D_00, "float", OPS_READ),   // t - 1
                     ops_arg_dat(dat_ut[(t + 1) % 3], 1, S2D_00, "float", OPS_WRITE),  // t + 1
                     ops_arg_gbl(&dt, 1, "float", OPS_READ), ops_arg_idx());

        for (int p_src = p_src_m; p_src <= p_src_M; p_src += 1) {
            float r0 = (int)(floor(-1.0e-1F * o_x + 1.0e-1F * src_coords[p_src]));
            ii_src[0] = r0;
            float r1 = (int)(floor(-1.0e-1F * o_y + 1.0e-1F * src_coords[p_src + 1]));
            ii_src[1] = r1;
            ii_src[2] = r1 + 1;
            ii_src[3] = r0 + 1;

            p[0] = (float)(-1.0e+1F * r0 - o_x + src_coords[p_src]);
            p[1] = (float)(-1.0e+1F * r1 - o_y + src_coords[p_src + 1]);
            // ops_printf("r0=%f r1=%f - o_x=%f o_y=%f - src_coords=%f %f\n", r0, r1, o_x, o_y, src_coords[p_src],
            //            src_coords[p_src + 1]);
            // ops_printf("p[0]=%f p[1]%f\n", p[0], p[1]);

            int source_range[] = {ii_src[1] + SPACE_ORDER, ii_src[2] + SPACE_ORDER + 1, ii_src[0] + SPACE_ORDER,
                                  ii_src[3] + SPACE_ORDER + 1};

            // ops_printf("src_range: %d %d %d %d\n", ii_src[1] + SPACE_ORDER, ii_src[2] + SPACE_ORDER + 1,
            //            ii_src[0] + SPACE_ORDER, ii_src[3] + SPACE_ORDER + 1);
            // printf("p_src=%d\tsrc_coords=%f\tsrc=%f\n", p_src, src_coords[p_src], src[t]);
            ops_par_loop(source_kernel, "source_kernel", grid, 2, source_range,
                         ops_arg_dat(dat_ut[(t + 1) % 3], 1, S2D_00, "float", OPS_WRITE),
                         ops_arg_dat(dat_m, 1, S2D_m2_m2, "float", OPS_READ),
                         ops_arg_gbl(&src[t], 1, "float", OPS_READ), ops_arg_gbl(p, 2, "float", OPS_READ),
                         ops_arg_gbl(ii_src, 4, "int", OPS_READ), ops_arg_gbl(x_limits, 2, "int", OPS_READ),
                         ops_arg_gbl(y_limits, 2, "int", OPS_READ), ops_arg_gbl(&dt, 1, "float", OPS_READ),
                         ops_arg_idx());
        }

        // if (t == 2) ops_print_dat_to_txtfile(dat_ut[(t + 1) % 3], "output-ops/u-partial-2.txt");
        // if (t == 10) ops_print_dat_to_txtfile(dat_ut[(t + 1) % 3], "output-ops/u-partial-10.txt");
        // if (t == 50) ops_print_dat_to_txtfile(dat_ut[(t + 1) % 3], "output-ops/u-partial-50.txt");
        // if (t == 100) ops_print_dat_to_txtfile(dat_ut[(t + 1) % 3], "output-ops/u-partial-100.txt");
        // ops_fetch_dat(dat_ut[(t + 1) % 3], u[(t + 1) % 3]);
        // array2d_to_file("output-ops/u.txt", u[(t + 1) % 3], size_u[0], size_u[1]);

        for (int p_rec = p_rec_m; p_rec <= p_rec_M; p_rec += 1) {
            // ops_printf("rec_coords=%f %f\n", rec_coords[p_rec], rec_coords[p_rec + (p_rec_M + 1)]);
            float r22 = (int)(floor(-1.0e-1F * o_x + 1.0e-1F * rec_coords[p_rec * 2]));
            ii_rec[0] = r22;
            float r23 = (int)(floor(-1.0e-1F * o_y + 1.0e-1F * rec_coords[p_rec * 2 + 1]));
            ii_rec[1] = r23;
            ii_rec[2] = r23 + 1;
            ii_rec[3] = r22 + 1;
            p[0] = (float)(-1.0e+1F * r22 - o_x + rec_coords[p_rec * 2]);
            p[1] = (float)(-1.0e+1F * r23 - o_y + rec_coords[p_rec * 2 + 1]);

            int receiver_range[] = {ii_rec[1] + SPACE_ORDER, ii_rec[2] + SPACE_ORDER + 1, ii_rec[0] + SPACE_ORDER,
                                    ii_rec[3] + SPACE_ORDER + 1};

            // ops_printf("ii_rec=[%d %d %d %d] p=[%f %f]\n", ii_rec[0], ii_rec[1], ii_rec[2], ii_rec[3], p[0], p[1]);

            ops_par_loop(receiver_kernel, "receiver_kernel", grid, 2, receiver_range,
                         ops_arg_reduce(red_sum, 1, "float", OPS_INC),
                         ops_arg_dat(dat_ut[t % 3], 1, S2D_00, "float", OPS_READ), ops_arg_gbl(p, 2, "float", OPS_READ),
                         ops_arg_gbl(ii_rec, 4, "int", OPS_READ), ops_arg_gbl(x_limits, 2, "int", OPS_READ),
                         ops_arg_gbl(y_limits, 2, "int", OPS_READ), ops_arg_idx());

            ops_reduction_result(red_sum, &rec[t * (p_rec_M + 1) + p_rec]);
            // ops_printf("time %d rec Position %d -> Value %e \n", t, t * (p_rec_M + 1) + p_rec, rec[t * (p_rec_M + 1)
            // + p_rec]);
        }
        // printf("\n");
    }

    ops_print_dat_to_txtfile(dat_ut[(T_INTERVALS + 1) % 3], "output-ops/u.txt");
    array2d_to_file("output-ops/rec.txt", rec, T_INTERVALS + 2, p_rec_M + 1);

    ops_exit();

    free(m);
    free(damp);
    free(src);
    free(rec);
    free(u);
    free(src_coords);
    free(rec_coords);

    return 0;
}