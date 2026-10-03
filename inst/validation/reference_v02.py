"""Reference values for deconflate 0.2 (hazard ratios, attributable risk,
method comparison). Run from this directory: python3 reference_v02.py
Requires numpy and scipy. Values printed here are hard-coded in
tests/testthat/test-hazard-ratio.R and test-compare.R.
"""
import itertools
from math import factorial
import numpy as np
from scipy.optimize import brentq, fsolve
from python_reference import or_to_joint, excess_matrix, ipf

np.set_printoptions(precision=10, suppress=True)


def hr_methods(P, OR, HR):
    E = excess_matrix(P, OR)
    A = np.eye(len(P)) + E.T
    m = HR - 1
    pub = 1 + np.where(m == 0, 0, m ** 2 / (m + (A - np.eye(len(P))) @ m))
    b1 = np.linalg.solve(A, np.log(HR))
    cells, p = ipf(P, OR, iters=4000)

    def crude(beta):
        w = p * np.exp(cells @ beta)
        n1 = (cells * w[:, None]).sum(0)
        n0 = w.sum() - n1
        return np.log(n1 / P) - np.log(n0 / (1 - P))

    b = fsolve(lambda x: crude(x) - np.log(HR), b1, xtol=1e-14)
    assert np.allclose(crude(b), np.log(HR), atol=1e-10)
    return pub, np.exp(b1), np.exp(b), cells, p


def attributable(cells, p, hr, r, shapley=True):
    beta = np.log(hr)
    rel = np.exp(cells @ beta)
    h0 = brentq(lambda h: (p * (1 - np.exp(-h * rel))).sum() - r, 1e-12, 50, xtol=1e-15)
    r0 = 1 - np.exp(-h0)
    out = {'disease_free': r0, 'attributable': r - r0}
    if shapley:
        n = cells.shape[1]
        phi = np.zeros(n)
        loss = lambda x: 1 - np.exp(-h0 * np.exp(x @ beta)) - r0
        for c, pr in zip(cells, p):
            present = np.flatnonzero(c)
            k = len(present)
            for j in present:
                others = [q for q in present if q != j]
                for size in range(k):
                    for S in itertools.combinations(others, size):
                        x = np.zeros(n); x[list(S)] = 1
                        xj = x.copy(); xj[j] = 1
                        w = factorial(size) * factorial(k - size - 1) / factorial(k)
                        phi[j] += pr * w * (loss(xj) - loss(x))
        out['shapley'] = phi
    return out


if __name__ == "__main__":
    # 1. Supplement example with culling hazard ratios 1.5, 2.0, 1.3
    P = np.array([.10, .15, .20])
    OR = np.ones((3, 3)); OR[0, 1] = OR[1, 0] = 2; OR[1, 2] = OR[2, 1] = 3
    HR = np.array([1.5, 2.0, 1.3])
    pub, sim, glob, cells, p = hr_methods(P, OR, HR)
    print('supplement HR published   ', repr(pub))
    print('supplement HR simultaneous', repr(sim))
    print('supplement HR global      ', repr(glob))
    a = attributable(cells, p, glob, 0.25)
    print('supplement attributable r=0.25:', repr(a['disease_free']), repr(a['attributable']),
          'shapley', repr(a['shapley']), 'sum', a['shapley'].sum())
    naive = (P * (1 - np.exp(-np.log(1 / (1 - 0.25)) * 1))).sum()

    # 2. Global dairy (analysis inputs), culling as hazard ratios
    D = ['CK', 'CM', 'DA', 'DYS', 'LAM', 'MET', 'MF', 'OC', 'PTB', 'RP', 'SCK', 'SCM']
    inc = dict(CK=0.0306287320599544, CM=0.305081708702075, DA=0.0216266798689032,
               DYS=0.0611107185570655, LAM=0.254288506511192, MET=0.0958013664473386,
               MF=0.0241829408262443, OC=0.112683551256031, PTB=0.103964073833149,
               RP=0.12346845541711, SCK=0.479052187708235, SCM=0.409414831921559)
    P = np.array([inc[d] if d in ('PTB', 'SCM') else 1 - np.exp(-inc[d]) for d in D])
    ors = {"CK:CM": 2.13, "CK:LAM": 1.65, "CK:MF": 1.60, "CK:OC": 1.97, "CK:RP": 1.55,
           "CK:SCK": 6.95, "CK:SCM": 2.40, "CM:LAM": 2.10, "CM:PTB": 1.89, "CM:RP": 2.70,
           "CM:SCK": 1.64, "CM:SCM": 3.05, "DA:CM": 3.45, "DA:MF": 2.50, "DA:RP": 3.50,
           "DA:SCK": 3.87, "DA:SCM": 3.60, "DYS:LAM": 2.09, "DYS:OC": 0.40, "DYS:RP": 2.74,
           "LAM:OC": 2.63, "LAM:PTB": 2.70, "LAM:RP": 1.50, "LAM:SCK": 2.01, "MET:CK": 2.42,
           "MET:CM": 2.30, "MET:DA": 3.40, "MET:DYS": 2.95, "MET:LAM": 6.10, "MET:MF": 1.50,
           "MET:OC": 1.94, "MET:RP": 3.53, "MET:SCK": 1.94, "MF:DYS": 9.70, "MF:LAM": 3.60,
           "MF:RP": 2.40, "RP:OC": 2.18, "SCK:RP": 1.52}
    OR = np.ones((12, 12))
    for k, v in ors.items():
        x, y = k.split(':'); OR[D.index(x), D.index(y)] = OR[D.index(y), D.index(x)] = v
    HR = np.array([1.5001, 2.3, 2.851179, 1.258143, 1.744976, 1.116444, 2.999886, 1.62,
                   2.310508, 1.599928, 1.92, 1.449996])
    pub, sim, glob, cells, p = hr_methods(P, OR, HR)
    print('\nglobal dairy HR published   ', repr(pub))
    print('global dairy HR simultaneous', repr(sim))
    print('global dairy HR global      ', repr(glob))
    for r in (0.2366,):
        a = attributable(cells, p, glob, r, shapley=False)
        print('global dairy attributable r=%.4f:' % r, repr(a['disease_free']), repr(a['attributable']))

    # 3. UK dairy 2022: total values per method (published, simultaneous)
    D = ['CO', 'DA', 'DYS', 'FAS', 'GIN', 'LAM', 'MAS', 'MET', 'MF', 'NEO', 'PTB', 'RP', 'SCK']
    Pd = dict(CO=.09, DA=.03, DYS=.02, FAS=.10, GIN=.21, LAM=.30, MAS=.30, MET=.10, MF=.08,
              NEO=.15, PTB=.07, RP=.05, SCK=.22)
    P = np.array([Pd[d] for d in D])
    ors = "RP:MET 6.20;DA:SCK 4.25;RP:DYS 4.10;MET:DA 3.40;MET:DYS 3.20;LAM:PTB 2.70;MF:DA 2.50;MAS:MET 2.30;DA:RP 2.20;MAS:DA 2.10;SCK:MF 2.10;LAM:SCK 2.01;MAS:MF 1.90;MAS:PTB 1.89;MAS:CO 1.65;MAS:SCK 1.64;SCK:CO 1.60;MET:SCK 1.40;SCK:RP 1.20"
    OR = np.ones((13, 13))
    for s in ors.split(';'):
        pr, v = s.split(); x, y = pr.split(':')
        OR[D.index(x), D.index(y)] = OR[D.index(y), D.index(x)] = float(v)
    Y = dict(CO=0, DA=4.04, DYS=4.05, FAS=7.33, GIN=3.28, LAM=5.54, MAS=4.57, MET=3.95, MF=0.41,
             NEO=4.20, PTB=5.90, RP=7.38, SCK=3.05)
    F = dict(CO=11.26, DA=0, DYS=6.96, FAS=0, GIN=1.20, LAM=12.47, MAS=0, MET=4.74, MF=0,
             NEO=7.21, PTB=5.79, RP=2.74, SCK=1.50)
    H = dict(DA=3.83, LAM=3.40, MAS=2.78, MF=2.50, PTB=2.40, MET=2.20, SCK=2.10, DYS=1.90,
             NEO=1.60)
    c = 0.27
    ex = np.array([(lambda p11: p11 / Pd[d] - (c - p11) / (1 - Pd[d]))(or_to_joint(c, Pd[d], H.get(d, 1)))
                   for d in D])
    E = excess_matrix(P, OR)
    A = np.eye(13) + E.T
    outs = [(np.array([Y[d] for d in D]) / 100, 'dec', 8737, 0.3022),
            (np.array([F[d] for d in D]) / 100, 'inc', 401, 13 * 0.3022),
            (ex, 'inc', 27, 1335.36 / 100)]
    for name in ('published', 'simultaneous'):
        tot = 71.09
        for m, dirn, x, price in outs:
            if name == 'published':
                adj = np.where(m == 0, 0, m ** 2 / np.where(m == 0, 1, m + (A - np.eye(13)) @ m))
            else:
                adj = np.linalg.solve(A, m)
            L = adj @ P
            xh = x / (1 - L) if dirn == 'dec' else x / (1 + L)
            tot += abs(xh - x) * price
        print('UK 2022 total value %-12s' % name, repr(tot))
