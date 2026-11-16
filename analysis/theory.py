"""Back-of-envelope theory for the run-based (patience / min-anchors) ANI method.

Model: sketched k-mers along genome A form a sequence of Bernoulli trials.
  - inside a homologous region: P(match) = q = ANI^k
  - outside (non-homologous):   P(match) = s  (spurious, ~ 2|B| / 4^k for random seq)
Gaps G between consecutive matches ~ Geometric(q) on {0,1,...}: P(G=g) = (1-q)^g q.
Patience x: a region closes after x consecutive misses.
"""
import math, random

C = 125
G_LEN = 5_000_000
ANIS = [0.99, 0.95, 0.90, 0.85, 0.80]
KS = [15, 16, 17, 19, 21]


def q_of(a, k):
    return a ** k


def s_of(k, g=G_LEN):
    return 1 - math.exp(-2 * g / 4 ** k)


def cond_gap_mean(q, x):
    """E[G | G < x] for G ~ Geom(q) on {0..}."""
    r = 1 - q
    num = r * (1 - x * r ** (x - 1) + (x - 1) * r ** x)
    return num / (q * (1 - r ** x))


def naive_ani(q, x, k):
    """ANI from pooled m/n inside regions (gaps truncated at x)."""
    qhat = 1 / (1 + cond_gap_mean(q, x))
    return qhat ** (1 / k)


def af_loss(q, x):
    """Fraction of sketched k-mers in a homologous region lost to false breaks."""
    r = 1 - q
    return r ** x * (q * x + r)


def min_x_for_loss(q, eps=0.01):
    x = 1
    while af_loss(q, x) > eps:
        x += 1
    return x


print("=== 1. match prob q = ANI^k  [mean bp between matches at c=%d] ===" % C)
print("ANI   " + "".join(f"     k={k:<10}" for k in KS))
for a in ANIS:
    row = ""
    for k in KS:
        q = q_of(a, k)
        row += f"  {q:6.4f} [{C / q / 1000:5.2f}kb]"
    print(f"{a:.2f} {row}")

print("\n=== 2. spurious match prob s per sketched k-mer (random 5 Mbp, both strands) ===")
for k in KS:
    print(f"k={k}: s={s_of(k):.2e}")

print("\n=== 3. patience x (sketched k-mers) so AF loss from false breaks <= 1% ===")
print("       (x*c bp = smallest non-homologous island that can be detected)")
print("ANI   " + "".join(f"   k={k:<12}" for k in KS))
for a in ANIS:
    row = ""
    for k in KS:
        q = q_of(a, k)
        x = min_x_for_loss(q)
        row += f"  x={x:4d} ({x * C / 1000:5.1f}kb)"
    print(f"{a:.2f} {row}")

print("\n=== 4. naive ANI bias (pooled m/n) when fixed x is too small, k=17 ===")
k = 17
for x in [10, 20, 40, 80, 160]:
    row = ""
    for a in ANIS:
        row += f"  {a:.2f}->{naive_ani(q_of(a, k), x, k):.4f}"
    print(f"x={x:4d}:{row}")

print("\n=== 5. spurious-region pressure s*x at x tuned for ANI_min (want << 0.1) ===")
for amin in [0.95, 0.90, 0.85, 0.80]:
    row = ""
    for k in KS:
        x = min_x_for_loss(q_of(amin, k))
        row += f"  k={k}:{s_of(k) * x:6.3f}"
    print(f"ANI_min={amin:.2f}:{row}")

print("\n=== 6. expected false regions/genome ~ n*s*(s*x)^(y-1), n=G/c, x tuned for ANI_min=0.85 ===")
n = G_LEN / C
for k in [15, 16, 17, 19]:
    x = min_x_for_loss(q_of(0.85, k))
    s = s_of(k)
    row = "".join(f"  y={y}:{n * s * min(1, s * x) ** (y - 1):9.3g}" for y in [2, 3, 4, 5])
    print(f"k={k} (x={x}):{row}")

print("\n=== 7. smallest contig (kb) expected to yield >= y anchors: y*c/q ===")
k = 17
for a in ANIS:
    q = q_of(a, k)
    print(f"ANI={a:.2f}: " + "  ".join(f"y={y}:{y * C / q / 1000:5.1f}kb" for y in [3, 4, 5]))

print("\n=== 8. adaptive patience x(q) = ceil(ln alpha / ln(1-q)), alpha=1e-3 ===")
for k in [15, 17, 19]:
    print(f"k={k}: " + "  ".join(
        f"{a:.2f}:{math.ceil(math.log(1e-3) / math.log(1 - q_of(a, k))):4d}" for a in ANIS))

print("\n=== 9. variance: per-bp Poisson rate vs per-sketched-kmer count, ratio Var_bp/Var_cnt = 1/(1-q) ===")
for a in [0.999, 0.99, 0.95, 0.90]:
    q = q_of(a, 15)
    print(f"ANI={a}: q={q:.3f}  variance ratio={1 / (1 - q):6.1f}x")


# --- Monte Carlo check of section 4 + truncated-geometric MLE fix --------------
def simulate(q, x, n=400_000, seed=1):
    rnd = random.Random(seed)
    bits = [rnd.random() < q for _ in range(n)]
    regions, cur, miss = [], None, 0
    for i, b in enumerate(bits):
        if b:
            if cur is None:
                cur = [i, i, 1]
            else:
                cur[1] = i
                cur[2] += 1
            miss = 0
        elif cur is not None:
            miss += 1
            if miss >= x:
                regions.append(cur)
                cur, miss = None, 0
    if cur:
        regions.append(cur)
    m = sum(r[2] for r in regions)
    span = sum(r[1] - r[0] + 1 for r in regions)
    gaps = span - m  # misses inside regions
    ngaps = m - len(regions)  # number of inter-match gaps
    return m, span, gaps, ngaps


def trunc_geom_mle(mean_gap, x):
    lo, hi = 1e-6, 1 - 1e-9
    for _ in range(100):
        mid = (lo + hi) / 2
        if cond_gap_mean(mid, x) > mean_gap:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


print("\n=== 10. Monte Carlo (k=17): true ANI vs naive m/n vs truncated-geometric MLE ===")
k = 17
for a, x in [(0.85, 20), (0.85, 40), (0.90, 15), (0.95, 6)]:
    q = q_of(a, k)
    m, span, gaps, ngaps = simulate(q, x)
    naive = (m / span) ** (1 / k)
    mle = trunc_geom_mle(gaps / ngaps, x) ** (1 / k)
    print(f"ANI={a} x={x:3d}: naive={naive:.4f} (theory {naive_ani(q, x, k):.4f})  trunc-MLE={mle:.4f}")
