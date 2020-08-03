#include <stdlib.h>
#include "math.h"
#include "source.h"

float *ricker_source(float fpeak, int length, float dt)
{
    float *source;
    int ti;
    float x, xx, t, tdelay = 1.0 / fpeak;

    // Allocate source array
    source = (float *)malloc(sizeof(float) * length);

    for (ti = 0; ti < length; ti++)
    {
        t = ti * dt;
        x = M_PI * fpeak * (t - tdelay);
        xx = x * x;
        source[ti] = exp(-xx) * (1.0 - 2.0 * xx);
    }

    return source;
}
