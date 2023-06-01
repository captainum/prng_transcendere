#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <gmp.h>
#include <math.h>


void func_w(mpz_t u[], mpz_t x[], mpz_t d, mpz_t b, mpz_t n, mpz_t res_up, mpz_t res_down, short fl)
{
    static mpz_t tmp0, tmp1, tmp2, tmp3, tmp4, dnx0, dnx1, dnx2, dnx3;
    static short tmp_inited = 0;
    if (!tmp_inited) {
        mpz_inits(tmp0, tmp1, tmp2, tmp3, tmp4, dnx0, dnx1, dnx2, dnx3, NULL);
        tmp_inited++;
    }

    // d * n + x[0]
    mpz_mul(tmp0, d, n);
    mpz_add(dnx0, tmp0, x[0]);

    // d * n + x[1]
    // mpz_mul(tmp0, d, n);
    mpz_add(dnx1, tmp0, x[1]);

    // d * n + x[2]
    // mpz_mul(tmp0, d, n);
    mpz_add(dnx2, tmp0, x[2]);

    // d * n + x[3]
    // mpz_mul(tmp0, d, n);
    mpz_add(dnx3, tmp0, x[3]);

    // u[0] * (d * n + x[1]) * (d * n + x[2]) * (d * n + x[3])
    mpz_mul(tmp1, u[0], dnx1);
    mpz_mul(tmp2, tmp1, dnx2);
    mpz_mul(tmp3, tmp2, dnx3);

    mpz_set(res_up, tmp3);
    //

    // u[1] * (d * n + x[0]) * (d * n + x[2]) * (d * n + x[3])
    mpz_mul(tmp1, u[1], dnx0);
    mpz_mul(tmp2, tmp1, dnx2);
    mpz_mul(tmp3, tmp2, dnx3);

    mpz_add(tmp4, res_up, tmp3);

    mpz_set(res_up, tmp4);
    //

    // u[2] * (d * n + x[0]) * (d * n + x[1]) * (d * n + x[3])
    mpz_mul(tmp1, u[2], dnx0);
    mpz_mul(tmp2, tmp1, dnx1);
    mpz_mul(tmp3, tmp2, dnx3);

    mpz_add(tmp4, res_up, tmp3);

    mpz_set(res_up, tmp4);
    //

    // u[3] * (d * n + x[0]) * (d * n + x[1]) * (d * n + x[2])
    mpz_mul(tmp1, u[3], dnx0);
    mpz_mul(tmp2, tmp1, dnx1);
    mpz_mul(tmp3, tmp2, dnx2);

    mpz_add(tmp4, res_up, tmp3);

    mpz_set(res_up, tmp4);
    //

    mpz_mul(tmp0, dnx0, dnx1);
    mpz_mul(tmp1, dnx2, dnx3);
    mpz_mul(tmp2, tmp0, tmp1);

    mpz_pow_ui(tmp3, b, mpz_get_ui(n));
    mpz_mul(res_down, tmp2, tmp3);

    if (fl) {
        mpz_clears(tmp0, tmp1, tmp2, tmp3, tmp4, dnx0, dnx1, dnx2, dnx3, NULL);
    }
}

void seq(mpz_t u[], mpz_t x[], mpz_t d, mpz_t bb, mpz_t b, uint64_t k, short **res) {
    static mpz_t nn, tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
    mpz_inits(nn, tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, NULL);

    mpz_t *w_up, *w_down;
    w_up = calloc(k + 1, sizeof(mpz_t));
    w_down = calloc(k + 1, sizeof(mpz_t));

    for (uint64_t i = 0; i < k + 1; i++) {
        mpz_inits(w_up[i], w_down[i], NULL);
    }

    uint64_t n = 0, n0 = 0;

    while (1) {
        if (n == k + 1) {
            w_up = realloc(w_up, (n + 1) * sizeof(mpz_t));
            w_down = realloc(w_down, (n + 1) * sizeof(mpz_t));

            mpz_inits(w_up[n], w_down[n], NULL);
        }

        mpz_init_set_ui(nn, n);
        func_w(u, x, d, bb, nn, w_up[n], w_down[n], 0);

        mpz_pow_ui(tmp0, b, n);
        mpz_mul(tmp1, w_up[n], tmp0);

        if (mpz_cmp(w_down[n], tmp1) > 0) {
            n0 = n;
            break;
        }
        else {
            n++;
        }
    }

    uint64_t h = n0;
    if (k + 1 > h) {
        h = k + 1;
        w_up = realloc(w_up, (h + 1) * sizeof(mpz_t));
        w_down = realloc(w_down, (h + 1) * sizeof(mpz_t));

        for (uint64_t i = n0 + 1; i < h + 1; i++) {
            mpz_inits(w_up[i], w_down[i], NULL);
        }
    }

    for (uint64_t i = n + 1; i < h + 1; i++) {
        mpz_init_set_ui(nn, i);
        func_w(u, x, d, bb, nn, w_up[i], w_down[i], 0);
    }

    mpz_t *delta_up, *delta_down, *alpha_up, *alpha_down;
    delta_up = calloc(h + 1, sizeof(mpz_t));
    delta_down = calloc(h + 1, sizeof(mpz_t));
    alpha_up = calloc(h + 1, sizeof(mpz_t));
    alpha_down = calloc(h + 1, sizeof(mpz_t));

    short *a;
    a = calloc(h + 1, sizeof(short));

    short c = 0;
    for (uint64_t n = 0; n < h + 1; n++) {
        if (c == 0) {
            {
                mpz_pow_ui(tmp0, b, n);
                mpz_mul(tmp1, tmp0, w_up[n]);
                mpz_init_set(alpha_up[n], tmp1);
            }

            {
                mpz_init_set(alpha_down[n], w_down[n]);
            }

            c = 1;
        }
        else {
            {
                mpz_mul(tmp1, b, delta_up[n-1]);
                mpz_mul(tmp2, tmp1, w_down[n]);

                mpz_pow_ui(tmp3, b, n);
                mpz_mul(tmp4, tmp3, delta_down[n-1]);
                mpz_mul(tmp5, tmp4, w_up[n]);

                mpz_add(tmp6, tmp2, tmp5);

                mpz_init_set(alpha_up[n], tmp6);
            }

            {
                mpz_mul(tmp0, delta_down[n-1], w_down[n]);

                mpz_init_set(alpha_down[n], tmp0);
            }
        }

        {
            mpz_t q;
            mpz_init(q);

            mpz_fdiv_q(q, alpha_up[n], alpha_down[n]);

            a[n] = mpz_get_ui(q);

            mpz_clear(q);
        }

        {
            mpz_mul_si(tmp0, alpha_down[n], a[n]);
            mpz_sub(tmp1, alpha_up[n], tmp0);
            mpz_init_set(delta_up[n], tmp1);
        }

        {
            mpz_init_set(delta_down[n], alpha_down[n]);
        }
    }

    uint64_t remain_size = 0;
    uint64_t old_sz = h + 1;

    while (mpz_cmp_ui(b, a[h] + 1) <= 0) {
        h++;

        if (!remain_size) {
            w_up = realloc(w_up, (old_sz * 2) * sizeof(mpz_t));
            w_down = realloc(w_down, (old_sz * 2) * sizeof(mpz_t));

            delta_up = realloc(delta_up, (old_sz * 2) * sizeof(mpz_t));
            delta_down = realloc(delta_down, (old_sz * 2) * sizeof(mpz_t));
            alpha_up = realloc(alpha_up, (old_sz * 2) * sizeof(mpz_t));
            alpha_down = realloc(alpha_down, (old_sz * 2) * sizeof(mpz_t));

            remain_size = (old_sz - 1) + 1;
            old_sz = old_sz * 2;
        }

        mpz_init_set_ui(nn, h);

        func_w(u, x, d, bb, nn, w_up[h], w_down[h], 0);

        {
            {
                mpz_mul(tmp1, b, delta_up[h-1]);
                mpz_mul(tmp1, tmp1, w_down[h]);

                mpz_pow_ui(tmp2, b, h);
                mpz_mul(tmp2, tmp2, delta_down[h-1]);
                mpz_mul(tmp2, tmp2, w_up[h]);
                mpz_add(tmp0, tmp1, tmp2);

                mpz_init_set(alpha_up[h], tmp0);
            }

            {
                mpz_mul(tmp0, delta_down[h-1], w_down[h]);

                mpz_init_set(alpha_down[h], tmp0);
            }
        }

        {
            {
                mpz_t q;
                mpz_init(q);

                mpz_fdiv_q(q, alpha_up[h], alpha_down[h]);

                a[h] = mpz_get_ui(q);

                mpz_clear(q);
            }

            {
                mpz_mul_si(tmp0, alpha_down[h], a[h]);
                mpz_sub(tmp1, alpha_up[h], tmp0);
                mpz_set(delta_up[h], tmp1);
            }

            {
                mpz_set(delta_down[h], alpha_down[h]);
            }
        }

        remain_size--;
    }

    for (int n = h - 1; n > 0; n--) {
        mpz_t q, r;
        mpz_inits(q, r, NULL);

        if ((a[n] < 0) || !mpz_cmp_si(b, a[n])) {
            mpz_set_si(tmp0, a[n]);

            mpz_fdiv_q(q, tmp0, b);
            mpz_fdiv_r(r, tmp0, b);

            a[n] = mpz_get_si(r);
            a[n-1] = a[n-1] + mpz_get_si(q);
        }

        mpz_clears(q, r, NULL);
    }

    for (unsigned int i = 0; i < h + 1; i++) {
        mpz_clears(w_up[i], w_down[i], delta_up[i], delta_down[i], alpha_up[i], alpha_down[i], NULL);
    }

    free(w_up);
    free(w_down);
    free(delta_up);
    free(delta_down);
    free(alpha_up);
    free(alpha_down);

    *res = a;
}

int main(int argc, char** argv) {
    int tmp;
    uint64_t tmp1;

    mpz_t u[4], x[4], d, bb, b;
    uint64_t k;
    short *a;

    FILE* f = fopen("input.txt", "r");

    if (!f) {
        exit(1);
    }

    for (int i = 0; i < 4; i++) {
        fscanf(f, "%d", &tmp);
        mpz_init_set_si(u[i], tmp);
    }

    for (int i = 0; i < 4; i++) {
        fscanf(f, "%d", &tmp);
        mpz_init_set_si(x[i], tmp);
    }

    fscanf(f, "%" PRIu64, &tmp1);
    mpz_init_set_ui(d, tmp1);

    fscanf(f, "%" PRIu64, &tmp1);
    mpz_init_set_ui(bb, tmp1);

    fscanf(f, "%" PRIu64, &tmp1);
    mpz_init_set_ui(b, tmp1);

    fscanf(f, "%" PRIu64, &k);

    fclose(f);
    
    seq(u, x, d, bb, b, k + 1, &a);

    f = fopen("output.txt", "w");

    if (!f) {
        exit(1);
    }

    for(uint64_t i = 1; i < k + 1; i++) {
        fprintf(f, "%d", a[i]);
    }

    fclose(f);

    for (int i = 0; i < 4; i++) {
        mpz_clears(u[i], x[i], NULL);
    }
    mpz_clears(d, bb, b, NULL);

    free(a);
    
    return 0;
}
