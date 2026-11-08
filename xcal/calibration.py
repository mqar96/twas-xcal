import numpy as np
from scipy.stats import norm, t as tdist
from scipy.optimize import brentq
from .core import standardize, panel_moments, jackknife_moments, satterthwaite_df


SIM_M_EFF = 312_000


class Panel:
    def __init__(self, G):
        self.Z, self.s, self.ok = standardize(G)
        self.n = self.Z.shape[0]
        self.R = self.Z.T @ self.Z / self.n

    def v(self, w):
        return np.where(self.ok, w * self.s, 0.0)

    def moments(self, w):
        return panel_moments(self.Z, self.v(w), self.ok)

    def jackknife(self, w, K=20, seed=0):
        return jackknife_moments(self.Z, self.s, np.where(self.ok, w, 0.0), K=K, seed=seed)


def true_var(v, R_gwas, c, intercept=1.0):
    Rv = R_gwas @ v
    return intercept * float(v @ Rv) + c * float(Rv @ Rv)


def type1(V, alpha):
    za = norm.isf(alpha / 2)
    return 2 * norm.sf(za / np.sqrt(np.asarray(V, float)))


def type1_t(V, df, alpha):
    thr = tdist.isf(alpha / 2, np.asarray(df, float))
    return 2 * norm.sf(thr / np.sqrt(np.asarray(V, float)))


def gc_lambda_analytic(V):
    V = np.asarray(V, float)
    f = lambda x: np.mean(2 * norm.cdf(np.sqrt(x / V)) - 1) - 0.5
    med = brentq(f, 1e-8, 1e4 * V.max())
    return med / 0.454936423119572


def type1_gc(V, alpha):
    lam = max(gc_lambda_analytic(V), 1.0)
    return type1(np.asarray(V, float) / lam, alpha), lam
