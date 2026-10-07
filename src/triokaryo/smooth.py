"""Smoothing of a track: the step fit by total-variation denoising, pooled windows of sites, and the null expectation of
the folded B-allele fraction.

The step fit solves  min_x  1/2 sum (y_i - x_i)^2 + lambda sum |x_{i+1} - x_i|  (the fused lasso's spatial part, the
taut string): a piecewise-constant reading of a noisy track in which every jump has to earn its height. Condat's direct
algorithm (IEEE Signal Processing Letters 20, 2013) finds the exact solution in one pass."""
import numpy as np

NA = float("nan")


def tv_denoise(y, lam):
    """The exact 1-D total-variation denoising of y (no NaN) with penalty lam; Condat's direct algorithm."""
    y = np.asarray(y, dtype=float)
    n = len(y)
    x = np.empty(n)
    if n == 0:
        return x
    if n == 1 or lam <= 0:
        x[:] = y
        return x
    k = k0 = kp = km = 0
    umin, umax = lam, -lam
    vmin, vmax = y[0] - lam, y[0] + lam
    while True:
        while k == n - 1:
            if umin < 0.0:                                      # a negative jump is needed at the end
                x[k0:km + 1] = vmin
                k0 = km + 1
                k = km = k0
                vmin = y[k]
                umin = lam
                umax = vmin + umin - vmax
            elif umax > 0.0:                                    # a positive one
                x[k0:kp + 1] = vmax
                k0 = kp + 1
                k = kp = k0
                vmax = y[k]
                umax = -lam
                umin = vmax + umax - vmin
            else:                                               # the last segment, at its mean
                vmin += umin / (k - k0 + 1)
                x[k0:] = vmin
                return x
        umin += y[k + 1] - vmin
        if umin < -lam:                                         # a negative jump is unavoidable: the segment ends at km
            x[k0:km + 1] = vmin
            k0 = km + 1
            k = kp = km = k0
            vmin = y[k]
            vmax = vmin + 2 * lam
            umin, umax = lam, -lam
            continue
        umax += y[k + 1] - vmax
        if umax > lam:                                          # a positive one: the segment ends at kp
            x[k0:kp + 1] = vmax
            k0 = kp + 1
            k = kp = km = k0
            vmax = y[k]
            vmin = vmax - 2 * lam
            umin, umax = lam, -lam
            continue
        k += 1                                                  # no jump yet: the running bounds move
        if umin >= lam:
            km = k
            vmin += (umin - lam) / (k - k0 + 1)
            umin = lam
        if umax <= -lam:
            kp = k
            vmax += (umax + lam) / (k - k0 + 1)
            umax = -lam


def tv_denoise_nan(y, lam):
    """The step fit over the finite values, NaN kept where y has none."""
    y = np.asarray(y, dtype=float)
    out = np.full(len(y), NA)
    ok = np.isfinite(y)
    if ok.sum():
        out[ok] = tv_denoise(y[ok], lam)
    return out


def tv_lambda(y, k=2.5):
    """A penalty for a track with noise sd from its first differences: k times that sd. A lone spike under 2 lambda is
    flattened; a plateau of m points keeps its height less 2 lambda / m."""
    v = np.asarray(y, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) < 5:
        return 0.0
    d = np.diff(v)
    sd = 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2)
    return float(k * sd) if np.isfinite(sd) and sd > 0 else 0.0


def windows_by_count(pos, w, max_gap=3_000_000):
    """Consecutive sites in windows of w, a window never spanning a gap of more than max_gap (a centromere), the last
    window of a run folded into the one before when it is under half full. Returns [(i0, i1)] over the sorted positions."""
    n = len(pos)
    if n == 0:
        return []
    cuts = [0] + [i + 1 for i in np.flatnonzero(np.diff(pos) > max_gap)] + [n]
    out = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        starts = list(range(a, b, w))
        for i, s in enumerate(starts):
            e = int(min(s + w, b))
            if i and e - s < w / 2:
                out[-1] = (out[-1][0], e)
            else:
                out.append((int(s), e))
    return out


_ABS_TABLE = {}


def folded_null(dp, max_table=400):
    """E|K/n - 1/2| for K ~ Binomial(n, 1/2), per site: what the folded B-allele fraction averages to at an ordinary
    heterozygous site of depth n (0.07 at 30x, 0.13 at 10x). Exact up to max_table, the normal approximation beyond."""
    dp = np.asarray(dp, dtype=int)
    out = np.full(len(dp), NA)
    if not _ABS_TABLE:
        t = np.zeros(max_table + 1)
        for n in range(1, max_table + 1):
            k = np.arange(n + 1)
            logp = np.cumsum(np.concatenate([[0.0], np.log(np.arange(n, 0, -1)) - np.log(np.arange(1, n + 1))])) - n * np.log(2.0)
            t[n] = float((np.abs(k / n - 0.5) * np.exp(logp)).sum())
        _ABS_TABLE["t"] = t
    t = _ABS_TABLE["t"]
    small = (dp >= 1) & (dp <= max_table)
    out[small] = t[dp[small]]
    big = dp > max_table
    out[big] = np.sqrt(2.0 / (np.pi * dp[big])) / 2.0
    return out
