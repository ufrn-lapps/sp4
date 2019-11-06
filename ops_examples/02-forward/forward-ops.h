// #ifndef FORWARD_OPS_H
// #define FORWARD_OPS_H

void propagation_kernel(const float *mx0, const float *dampx0, const float *utime0, const float *utime2, float *utime1,
                        const float *dt, const int *idx) {
    utime1[OPS_ACC4(0, 0)] = 2.0F *
                             (2.5e-1F * *dt * dampx0[OPS_ACC1(-3, -3)] * utime2[OPS_ACC3(0, 0)] +
                              (*dt * *dt) * (6.66666652e-3F * (utime0[OPS_ACC2(-1, 0)] + utime0[OPS_ACC2(0, -1)] +
                                                               utime0[OPS_ACC2(0, 1)] + utime0[OPS_ACC2(1, 0)]) -
                                             4.16666657e-4F * (utime0[OPS_ACC2(-2, 0)] + utime0[OPS_ACC2(0, -2)] +
                                                               utime0[OPS_ACC2(0, 2)] + utime0[OPS_ACC2(2, 0)]) -
                                             2.49999994e-2F * utime0[OPS_ACC2(0, 0)]) -
                              5.0e-1F * mx0[OPS_ACC0(-2, -2)] * utime2[OPS_ACC3(0, 0)] +
                              1.0F * mx0[OPS_ACC0(-2, -2)] * utime0[OPS_ACC2(0, 0)]) /
                             (5.0e-1F * *dt * dampx0[OPS_ACC1(-3, -3)] + 1.0F * mx0[OPS_ACC0(-2, -2)]);

    // if ((idx[0] == 90 && idx[1] == 80) || (idx[0] == 88 && idx[1] == 80))
    // if (idx[1] == 6) {
    //     ops_printf("damp=%f\tu[%d][%d]=%e\n", dampx0[OPS_ACC1(-3, -3)], idx[1], idx[0], utime1[OPS_ACC4(0, 0)]);
    // }
    // ops_printf("%f\n", utime1[OPS_ACC4(0, 0)]);
}

void source_kernel(float *u, const float *m, const float *src_value, const float *p, const int *ii_src,
                   const int *x_limits, const int *y_limits, const float *dt, const int *idx) {
    // ops_printf("src=%f\n", *src_value);
    if (ii_src[0] >= x_limits[0] - 1 && ii_src[1] >= y_limits[0] - 1 && ii_src[0] <= x_limits[1] + 1 &&
        ii_src[1] <= y_limits[1] + 1 && ii_src[0] + 4 == idx[1] && ii_src[1] + 4 == idx[0]) {
        float r6 = (*dt * *dt) * (1.0e-2F * p[0] * p[1] - 1.0e-1F * p[0] - 1.0e-1F * p[1] + 1) * *src_value /
                   m[OPS_ACC1(-2, -2)];
        u[OPS_ACC0(0, 0)] += r6;
        // ops_printf("dt=%f px=%f py=%f src_value=%f m=%f\n", *dt, p[1], p[0], *src_value, m[OPS_ACC1(-2, -2)]);
        // ops_printf("u[%d][%d] injected=%f\n", idx[1], idx[0], r6);
    }
    if (ii_src[0] >= x_limits[0] - 1 && ii_src[2] >= y_limits[0] - 1 && ii_src[0] <= x_limits[1] + 1 &&
        ii_src[2] <= y_limits[1] + 1 && ii_src[0] + 4 == idx[1] && ii_src[2] + 4 == idx[0]) {
        float r11 = (*dt * *dt) * (-1.0e-2F * p[0] * p[1] + 1.0e-1F * p[1]) * *src_value / m[OPS_ACC1(-2, -2)];
        u[OPS_ACC0(0, 0)] += r11;
        // ops_printf("u[%d][%d] injected=%f\n", idx[1], idx[0], r11);
    }
    if (ii_src[1] >= y_limits[0] - 1 && ii_src[3] >= x_limits[0] - 1 && ii_src[1] <= y_limits[1] + 1 &&
        ii_src[3] <= x_limits[1] + 1 && ii_src[3] + 4 == idx[1] && ii_src[1] + 4 == idx[0]) {
        float r16 = (*dt * *dt) * (-1.0e-2F * p[0] * p[1] + 1.0e-1F * p[0]) * *src_value / m[OPS_ACC1(-2, -2)];
        u[OPS_ACC0(0, 0)] += r16;
        // ops_printf("u[%d][%d] injected=%f\n", idx[1], idx[0], r16);
    }
    if (ii_src[2] >= y_limits[0] - 1 && ii_src[3] >= x_limits[0] - 1 && ii_src[2] <= y_limits[1] + 1 &&
        ii_src[3] <= x_limits[1] + 1 && ii_src[3] + 4 == idx[1] && ii_src[2] + 4 == idx[0]) {
        float r21 = 1.0e-2F * p[0] * p[1] * (*dt * *dt) * *src_value / m[OPS_ACC1(-2, -2)];
        u[OPS_ACC0(0, 0)] += r21;
        // ops_printf("u[%d][%d] in jected=%f\n", idx[1], idx[0], r21);
    }
}

void receiver_kernel(float *sum, const float *u, const float *p, const int *ii_rec, const int *x_limits,
                     const int *y_limits, const int *idx) {
    // ops_printf("%d %d %d %d %d %d\n", x_limits[0] - 1, y_limits[0] - 1, x_limits[1] + 1, y_limits[1] + 1, ii_rec[0] +
    // 4,
    //            ii_rec[1] + 4);
    if (idx[0] >= x_limits[0] - 1 && idx[1] >= y_limits[0] - 1 && idx[0] <= x_limits[1] + 1 &&
        idx[1] <= y_limits[1] + 1 && idx[1] == ii_rec[0] + 4 && idx[0] == ii_rec[1] + 4) {
        *sum += (1.0e-2F * p[0] * p[1] - 1.0e-1F * p[0] - 1.0e-1F * p[1] + 1) * u[OPS_ACC1(0, 0)];
        // ops_printf("u = %e p[0] = %f  p[1] = %f\n", u[OPS_ACC1(0, 0)], p[0], p[1]);
        // ops_printf("idx[0]=%d idx[1]=%d sum=%e\n", idx[0], idx[1], *sum);
    }
    if (idx[0] >= x_limits[0] - 1 && idx[1] >= y_limits[0] - 1 && idx[0] <= x_limits[1] + 1 &&
        idx[1] <= y_limits[1] + 1 && idx[1] == ii_rec[0] + 4 && idx[0] == ii_rec[2] + 4) {
        *sum += (-1.0e-2F * p[0] * p[1] + 1.0e-1F * p[1]) * u[OPS_ACC1(0, 0)];
        // ops_printf("idx[0]=%d idx[1]=%d sum=%e\n", idx[0], idx[1], *sum);
    }
    if (idx[0] >= y_limits[0] - 1 && idx[1] >= x_limits[0] - 1 && idx[0] <= y_limits[1] + 1 &&
        idx[1] <= x_limits[1] + 1 && idx[1] == ii_rec[3] + 4 && idx[0] == ii_rec[1] + 4) {
        *sum += (-1.0e-2F * p[0] * p[1] + 1.0e-1F * p[0]) * u[OPS_ACC1(0, 0)];
        // ops_printf("idx[0]=%d idx[1]=%d sum=%e\n", idx[0], idx[1], *sum);
    }
    if (idx[0] >= y_limits[0] - 1 && idx[1] >= x_limits[0] - 1 && idx[0] <= y_limits[1] + 1 &&
        idx[1] <= x_limits[1] + 1 && idx[1] == ii_rec[3] + 4 && idx[0] == ii_rec[2] + 4) {
        *sum += 1.0e-2F * p[0] * p[1] * u[OPS_ACC1(0, 0)];
        // ops_printf("idx[0]=%d idx[1]=%d sum=%e\n", idx[0], idx[1], *sum);
    }
}

// #endif