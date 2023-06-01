import math

def seq16(u: list, x: list, d: int, k: int):
    def func_w16(n: int):
        return (
            (
                u[0] * (d * n + x[1]) * (d * n + x[2]) * (d * n + x[3]) +
                u[1] * (d * n + x[0]) * (d * n + x[2]) * (d * n + x[3]) +
                u[2] * (d * n + x[0]) * (d * n + x[1]) * (d * n + x[3]) +
                u[3] * (d * n + x[0]) * (d * n + x[1]) * (d * n + x[2])
            ),
            (d * n + x[0]) * (d * n + x[1]) * (d * n + x[2]) * (d * n + x[3])  # * (16 ** n)
        )

    w_up, w_down = [], []
    n = 0

    while True:
        wu, wd = func_w16(n)

        w_up.append(wu)
        w_down.append(wd)

        if w_up[n] < w_down[n]:
            n0 = n
            break
        else:
            n += 1

    h = max(k + 1, n0)

    for i in range(n + 1, h + 1):
        wu, wd = func_w16(i)

        w_up.append(wu)
        w_down.append(wd)

    delta_up, delta_down = [], []
    alpha_up, alpha_down = [], []

    a = []

    for n in range(0, h + 1):
        if len(delta_up) == 0:
            alpha_up.append(16 * 0 * w_down[n] + 1 * w_up[n])
            alpha_down.append(1 * w_down[n])
        else:
            alpha_up.append(16 * delta_up[n - 1] * w_down[n] + delta_down[n - 1] * w_up[n])
            alpha_down.append(delta_down[n - 1] * w_down[n])

        a.append(alpha_up[n] // alpha_down[n])
        delta_up.append(alpha_up[n] - a[n] * alpha_down[n])
        delta_down.append(alpha_down[n])

    while a[h] >= 16 - 1:
        h += 1

        wu, wd = func_w16(n)

        w_up.append(wu)
        w_down.append(wd)

        alpha_up.append(16 * delta_up[h - 1] * w_down[h] + delta_down[h - 1] * w_up[h])
        alpha_down.append(delta_down[h - 1] * w_down[h])

        a.append(math.floor(alpha_up[h] / alpha_down[h]))
        delta_up.append(alpha_up[h] - a[h] * alpha_down[h])
        delta_down.append(alpha_down[h])

    for n in range(h - 1, 0, -1):
        if a[n] < 0 or a[n] >= 16:
            q = a[n] // 16
            r = pow(a[n], 1, 16)

            a[n] = r
            a[n - 1] = a[n - 1] + q

    return a


with open("input.txt") as f:
    u = [int(i) for i in next(f).split()]
    x = [int(i) for i in next(f).split()]
    d = int(next(f))
    bb = int(next(f))

    b = int(next(f))
    k = int(next(f)) + 1

a16 = seq16(u, x, d, k)
with open("tests/16/output_davasilchenko_16.txt", "w") as f:
    f.write(''.join(map(str, a16[1:k])))
