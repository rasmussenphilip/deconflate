# Reference values for the 2024 global dairy analysis inputs, as used by the
# published analysis code (Rasmussen et al. 2024, 1st revision):
#   * fixed de-conflation probabilities: 1 - exp(-mean incidence) for all
#     diseases except PTB (prevalence) and SCM (entered unconverted, 0.4094);
#   * odds ratios from Table 3 (PERT central value = mode; normals truncated
#     at 0 here, as in deconflate; the original code left them untruncated);
#   * unrounded impact distributions from the analysis input tables;
#   * culling entered as (HR - 1), de-conflated, then 1 added back.
# Mirrors example_global_dairy() and sampler_global_dairy() (inputs = "analysis").
# Run from this directory: python3 reference_2024_analysis.py [n_draws]
import sys
import numpy as np
from scipy.stats import norm

D = ['CK', 'CM', 'DA', 'DYS', 'LAM', 'MET', 'MF', 'OC', 'PTB', 'RP', 'SCK', 'SCM']
n = len(D)
inc = {'CK': 0.0306287320599544, 'CM': 0.305081708702075, 'DA': 0.0216266798689032,
       'DYS': 0.0611107185570655, 'LAM': 0.254288506511192, 'MET': 0.0958013664473386,
       'MF': 0.0241829408262443, 'OC': 0.112683551256031, 'PTB': 0.103964073833149,
       'RP': 0.12346845541711, 'SCK': 0.479052187708235, 'SCM': 0.409414831921559}
P = np.array([inc[d] if d in ('PTB', 'SCM') else 1 - np.exp(-inc[d]) for d in D])

# (d1, d2, type, params): f fixed; n normal(mean, sd) truncated at 0; p PERT(min, mode, max)
ORS = """CK CM p 1.2 2.13 3.4;CK LAM p 1.2 1.65 2.4;CK MF n 1.6 0.13;CK OC p 1.3 1.97 4.1;
CK RP p 1 1.55 1.9;CK SCK n 6.95 1.28;CK SCM n 2.4 0.41;CM LAM f 2.1;CM PTB n 1.89 0.2;
CM RP n 2.7 0.33;CM SCK n 1.64 0.2;CM SCM p 1.3 3.05 6.5;DA CM p 1.4 3.45 4.8;DA MF n 2.5 0.48;
DA RP p 1.6 3.5 4.6;DA SCK n 3.87 0.34;DA SCM n 3.6 1.35;DYS LAM n 2.09 0.26;DYS OC f 0.4;
DYS RP p 1.25 2.74 5.96;LAM OC n 2.63 1.44;LAM PTB n 2.7 1.22;LAM RP n 1.5 0.31;
LAM SCK n 2.01 0.2;MET CK p 1.2 2.42 10.4;MET CM p 1.2 2.3 3.8;MET DA p 1.6 3.4 7.6;
MET DYS p 0.98 2.95 9.72;MET LAM n 6.1 1.45;MET MF n 1.5 0.15;MET OC p 1.2 1.94 3;
MET RP p 1.8 3.53 6.52;MET SCK n 1.94 0.09;MF DYS n 9.7 1.3;MF LAM f 3.6;MF RP n 2.4 0.2;
RP OC p 1.78 2.18 2.57;SCK RP n 1.52 0.19"""

YIELD = {'CK': 'p 0.2351676 0.4321944 1.043482', 'CM': 'n 3.2499 0.7584468',
         'DA': 'p -1.451884 2.83693 9.190728', 'DYS': 'n 4.919088 0.9688162',
         'LAM': 'n 4.8061 0.8651519', 'MET': 'n 5.613085 1.350757',
         'MF': 'f 0.5365130859025376', 'OC': 'p 1.711971 3.747839 4.326951',
         'PTB': 'n 4.3 0.6683673', 'RP': 'n 4.198664 1.154555',
         'SCK': 'n 8.396472 1.185384', 'SCM': 'n 6.293184 1.200231'}
FERT = {'CK': 'n 1.445122 0.3624052', 'CM': 'n 8.42 2.424912', 'DA': 'n 1.082641 2.036204',
        'DYS': 'n 2.399177 0.9317402', 'LAM': 'p 1.190476 3.304898 10.71429',
        'MET': 'n 14.67308 8.535001', 'MF': 'p 2.032968 2.414949 3.095238',
        'OC': 'p 5.043478 9.685465 21.42857', 'PTB': 'n 5.349124 2.527459',
        'RP': 'n 6.760971 1.558483', 'SCK': 'n 1.122037 1.822338',
        'SCM': 'p -0.1242236 0.2636645 5.681529'}
CULL = {'CK': 'n 0.5001 0.10000976530231317', 'CM': 'n 1.3 0.17342740434782608',
        'DA': 'p 0 1.851179 6.9', 'DYS': 'p -0.4 0.258143 1.1',
        'LAM': 'n 0.744976 0.07449644729921788', 'MET': 'p -0.4 0.116444 0.5',
        'MF': 'n 1.999886 0.6012380424926813', 'OC': 'n 0.62 0.1591395716049383',
        'PTB': 'n 1.310508 0.2011405862857865', 'RP': 'n 0.599928 0.11954628295273287',
        'SCK': 'n 0.92 0.0855654625', 'SCM': 'n 0.449996 0.07758647609400302'}


def central(spec):
    t, *v = spec.split()
    v = list(map(float, v))
    return v[1] if t == 'p' else v[0]


def draw(spec, N, rng, lower=-np.inf):
    t, *v = spec.split()
    v = list(map(float, v))
    if t == 'f':
        return np.full(N, v[0])
    if t == 'n':
        Fl = norm.cdf(lower, v[0], v[1])
        return norm.ppf(Fl + rng.random(N) * (1 - Fl), v[0], v[1])
    mn, mode, mx = v
    a = 1 + 4 * (mode - mn) / (mx - mn)
    b = 1 + 4 * (mx - mode) / (mx - mn)
    return mn + (mx - mn) * rng.beta(a, b, N)


def excess(OR):
    """E[..., k, i] = P(k|i) - P(k|not i) for arrays of odds-ratio matrices."""
    N = OR.shape[0]
    E = np.zeros((N, n, n))
    for i in range(n):
        for k in range(n):
            if i == k:
                continue
            o = OR[:, i, k]
            pi, pk = P[i], P[k]
            a = (o - 1) * pk
            b = -((o - 1) * (pi + pk) + 1)
            c = o * pi
            with np.errstate(invalid='ignore', divide='ignore'):
                g = (-b - np.sqrt(b * b - 4 * a * c)) / (2 * a)
            p11 = np.where(np.abs(o - 1) < 1e-12, pi * pk, g * pk)
            E[:, k, i] = p11 / pi - (pk - p11) / (1 - pi)
    return E


def adjust(E, M):
    conf = np.einsum('ski,sk->si', E, M)
    with np.errstate(invalid='ignore', divide='ignore'):
        pub = np.where(M == 0, 0.0, M ** 2 / (M + conf))
    A = np.eye(n)[None] + np.transpose(E, (0, 2, 1))
    sim = np.linalg.solve(A, M[..., None])[..., 0]
    return pub, sim


pairs = []
for s in ORS.replace('\n', '').split(';'):
    d1, d2, *spec = s.split()
    pairs.append((D.index(d1), D.index(d2), ' '.join(spec)))
OUT = [('yield', YIELD, 100), ('fertility', FERT, 100), ('culling', CULL, 1)]

# 1. Central values (example_global_dairy(inputs = "analysis")).
OR0 = np.ones((1, n, n))
for i, k, spec in pairs:
    OR0[0, i, k] = OR0[0, k, i] = central(spec)
E0 = excess(OR0)
print('Probabilities:', ' '.join(f'{d}={p:.6f}' for d, p in zip(D, P)))
for lab, spec, div in OUT:
    M = np.array([[central(spec[d]) / div for d in D]])
    pub, sim = adjust(E0, M)
    print(f'central {lab} published  :', ', '.join(f'{x:.6f}' for x in pub[0]))
    print(f'central {lab} simultaneous:', ', '.join(f'{x:.6f}' for x in sim[0]))

# 2. Monte Carlo means (sampler_global_dairy(inputs = "analysis")).
N = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
rng = np.random.default_rng(2024)
OR = np.ones((N, n, n))
for i, k, spec in pairs:
    v = draw(spec, N, rng, lower=0)
    OR[:, i, k] = OR[:, k, i] = v
E = excess(OR)
T5 = {'yield': [.03, 1.36, 1.18, 3.48, 2.62, 2.87, .07, 2.59, 3.37, 2.30, 7.11, 5.58],
      'fertility': [.34, 6.09, .78, 1.11, 1.86, 11.22, 1.06, 9.03, 4.23, 3.74, .39, .04],
      'culling': [1.18, 1.90, 2.75, 1.18, 1.40, 1.03, 2.64, 1.51, 2.07, 1.29, 1.67, 1.25]}
for lab, spec, div in OUT:
    M = np.column_stack([draw(spec[d], N, rng) / div for d in D])
    pub, sim = adjust(E, M)
    scale, add = (100, 0) if div == 100 else (1, 1)
    mu = pub.mean(0) * scale + add
    se = pub.std(0) / np.sqrt(N) * scale
    print(f'\nMC {lab} (N={N}): disease  Table5  published_mean (MC SE)  median')
    for j, d in enumerate(D):
        print(f'  {d:4s} {T5[lab][j]:6.2f} {mu[j]:8.3f} ({se[j]:.3f}) {np.median(pub[:, j]) * scale + add:8.3f}')
