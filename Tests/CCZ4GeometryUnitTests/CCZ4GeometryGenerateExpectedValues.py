#!/usr/bin/env python3
# GRTeclyn
# Copyright 2022 The GRTL collaboration.
# Please refer to LICENSE in GRTeclyn's root directory.
"""Generate CCZ4GeometryExpectedValues.hpp for the CCZ4 geometry unit test.

This is a pure-Python replacement for the old Mathematica notebook
(Mathematica_comparison.nb). It draws a random set of CCZ4 variables and their
first/second derivatives, computes the geometric quantities that the unit test
checks (inverse metric, conformal Christoffel symbols and their contraction,
and the Z4 Ricci tensor + scalar), and writes both the input values and the
expected results into a single header.

The generated header contains two preprocessor-guarded sections so that the one
file can be #included at the two different points in the unit test:

    #define CCZ4_GEOMETRY_INPUT_VALUES
    #include "CCZ4GeometryExpectedValues.hpp"
    #undef  CCZ4_GEOMETRY_INPUT_VALUES

    #define CCZ4_GEOMETRY_EXPECTED_VALUES
    #include "CCZ4GeometryExpectedValues.hpp"
    #undef  CCZ4_GEOMETRY_EXPECTED_VALUES

Only the Python standard library is required (no numpy).

Index / storage conventions (matching the C++ Tensor classes):
  * The last index of a derivative tensor is always the derivative direction.
  * Symmetric index pairs are packed as
        (0,0)->0  (0,1)->1  (0,2)->2  (1,1)->3  (1,2)->4  (2,2)->5
    i.e. VAR_IDX0(i, j) = i + j + (1 if i*j != 0 else 0).
  * d1_h  is Sym12Rank3:      d1_h(pair(i,j), k)          = d_k h_{ij}
  * d1_Gamma is Rank2:        d1_Gamma(a, b)              = d_b Gamma^a
  * d2_chi is Sym12Rank2:     d2_chi(pair(i,j))           = d_i d_j chi
  * d2_h  is Sym12Sym34Rank4: d2_h(pair(i,j), pair(k,l))  = d_k d_l h_{ij}
"""

import argparse
import os
import random

DIM = 3

# Symmetric-pair index mapping, in packed order 0..5.
SYM_PAIRS = [(i, j) for i in range(DIM) for j in range(i, DIM)]


def var_idx0(i, j):
    """Packed index of the symmetric pair (i, j)."""
    return i + j + (1 if i * j != 0 else 0)


# --------------------------------------------------------------------------- #
# Small pure-Python linear algebra helpers
# --------------------------------------------------------------------------- #
def zeros(*shape):
    if len(shape) == 1:
        return [0.0] * shape[0]
    return [zeros(*shape[1:]) for _ in range(shape[0])]


def det3(m):
    """Determinant of a 3x3 matrix."""
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def inv3(m):
    """Inverse of a 3x3 matrix via the adjugate / determinant."""
    det = det3(m)
    cof = [[m[1][1] * m[2][2] - m[1][2] * m[2][1],
            m[0][2] * m[2][1] - m[0][1] * m[2][2],
            m[0][1] * m[1][2] - m[0][2] * m[1][1]],
           [m[1][2] * m[2][0] - m[1][0] * m[2][2],
            m[0][0] * m[2][2] - m[0][2] * m[2][0],
            m[0][2] * m[1][0] - m[0][0] * m[1][2]],
           [m[1][0] * m[2][1] - m[1][1] * m[2][0],
            m[0][1] * m[2][0] - m[0][0] * m[2][1],
            m[0][0] * m[1][1] - m[0][1] * m[1][0]]]
    return [[cof[i][j] / det for j in range(DIM)] for i in range(DIM)]


# --------------------------------------------------------------------------- #
# Random input generation
# --------------------------------------------------------------------------- #
def _random_sym():
    """Random symmetric 3x3 matrix with entries in [0, 1)."""
    m = [[random.random() for _ in range(DIM)] for _ in range(DIM)]
    for i in range(DIM):
        for j in range(i + 1, DIM):
            m[j][i] = m[i][j]
    return m


def unit_det_metric():
    """Random symmetric conformal metric with unit determinant, together with
    its consistent first and second derivatives (h, dh, d2h).

    The conformal metric of BSSN/CCZ4 is symmetric AND has det = 1, and the
    contracted-Christoffel derivative (needed for the general-form Ricci Z) is
    only well defined when dh, d2h are genuine derivatives of such a metric.

    A smooth symmetric field is built from random Taylor coefficients about the
    origin,
        g_ij(x) = S_ij + A[l]_ij x_l + 1/2 H[l][m]_ij x_l x_m,
    (S, A, H symmetric in i,j; H symmetric in l,m; S is diagonally boosted so g
    is positive definite) and normalised to unit determinant, h = phi g with
    phi = det(g)^(-1/3). Evaluating at x = 0 gives det(h) = 1 with dh, d2h the
    true derivatives of the det = 1 field:
        d_l phi     = -1/3 phi T_l,           T_l = tr(g^{-1} d_l g)
        d_l h_ij    = phi (d_l g_ij - 1/3 T_l g_ij)
        d_l d_m h_ij= phi [ (1/9 T_l T_m - 1/3 d_m T_l) g_ij
                            - 1/3 T_l d_m g_ij - 1/3 T_m d_l g_ij + d_l d_m g_ij ]
    """
    g = _random_sym()
    for i in range(DIM):
        g[i][i] += 2.0                          # ensure positive definiteness
    dg = [_random_sym() for _ in range(DIM)]    # dg[l]_ij = d_l g_ij
    d2g = [[None] * DIM for _ in range(DIM)]    # d2g[l][m]_ij, symmetric in l,m
    for l in range(DIM):
        for m in range(l, DIM):
            block = _random_sym()
            d2g[l][m] = block
            d2g[m][l] = block

    gu = inv3(g)
    phi = det3(g) ** (-1.0 / 3.0)

    # T_l = tr(g^{-1} d_l g) = d_l ln det g;   dT[m][l] = d_m T_l (sym in l,m)
    T = [sum(gu[a][b] * dg[l][a][b] for a in range(DIM) for b in range(DIM))
         for l in range(DIM)]
    dT = zeros(DIM, DIM)
    for m in range(DIM):
        for l in range(DIM):
            s1 = 0.0
            for a in range(DIM):
                for b in range(DIM):
                    for c in range(DIM):
                        for d in range(DIM):
                            s1 += gu[a][b] * dg[m][b][c] * gu[c][d] * dg[l][d][a]
            s2 = sum(gu[a][b] * d2g[m][l][a][b]
                     for a in range(DIM) for b in range(DIM))
            dT[m][l] = -s1 + s2

    h = [[phi * g[i][j] for j in range(DIM)] for i in range(DIM)]
    dh = zeros(DIM, DIM, DIM)
    for l in range(DIM):
        for i in range(DIM):
            for j in range(DIM):
                dh[l][i][j] = phi * (dg[l][i][j] - T[l] * g[i][j] / 3.0)
    d2h = zeros(DIM, DIM, DIM, DIM)
    for l in range(DIM):
        for m in range(DIM):
            coeff_g = T[l] * T[m] / 9.0 - dT[m][l] / 3.0
            for i in range(DIM):
                for j in range(DIM):
                    d2h[l][m][i][j] = phi * (
                        coeff_g * g[i][j]
                        - T[l] * dg[m][i][j] / 3.0
                        - T[m] * dg[l][i][j] / 3.0
                        + d2g[l][m][i][j])
    return h, dh, d2h


def generate_inputs():
    """Return a dict of the generated CCZ4 variables and derivatives.

    The conformal metric sector (h, d1_h, d2_h) is a consistent unit-determinant
    metric (see unit_det_metric). Every other quantity is an independent random
    field, exactly as the Mathematica notebook generated them -- in particular
    Gamma is the independent evolved connection variable (distinct from the
    metric-derived Gammabar), which is what makes the constraint vector
    Z_i = 1/2 (Gamma_i - Gammabar_i) nonzero.
    """
    chi = random.random()
    h, dh, d2h = unit_det_metric()

    gamma = [random.random() for _ in range(DIM)]         # Gamma^i
    dchi = [random.random() for _ in range(DIM)]          # d_i chi

    # dgamma[i][j] = d_i Gamma^j (no symmetry).
    dgamma = [[random.random() for _ in range(DIM)] for _ in range(DIM)]

    # d2chi[i][j] = d_i d_j chi, symmetric.
    d2chi = [[random.random() for _ in range(DIM)] for _ in range(DIM)]
    for i in range(DIM):
        for j in range(i + 1, DIM):
            d2chi[j][i] = d2chi[i][j]

    z_over_chi = [random.random() for _ in range(DIM)]

    return dict(chi=chi, h=h, gamma=gamma, dchi=dchi, dh=dh, dgamma=dgamma,
                d2chi=d2chi, d2h=d2h, z_over_chi=z_over_chi)


# --------------------------------------------------------------------------- #
# The geometry computation (port of the Mathematica notebook)
# --------------------------------------------------------------------------- #
def compute_expected(v):
    chi = v["chi"]
    h = v["h"]
    gamma = v["gamma"]
    dchi = v["dchi"]
    dh = v["dh"]
    dgamma = v["dgamma"]
    d2chi = v["d2chi"]
    d2h = v["d2h"]
    z_over_chi = v["z_over_chi"]

    dims = 4  # spacetime dimension; spatial factors use dims - 1 = 3
    hu = inv3(h)

    # Conformal Christoffel symbols: chris[m][b][c] = Gamma^m_{bc}
    #   = 1/2 h^{ma} (d_b h_{ca} + d_c h_{ba} - d_a h_{bc})
    chris = zeros(DIM, DIM, DIM)
    for m in range(DIM):
        for b in range(DIM):
            for c in range(DIM):
                s = 0.0
                for a in range(DIM):
                    s += hu[m][a] * (dh[b][c][a] + dh[c][b][a] - dh[a][b][c])
                chris[m][b][c] = 0.5 * s

    # Contracted Christoffel: chrisContr[c] = h^{ab} Gamma^c_{ab}
    chris_contr = zeros(DIM)
    for c in range(DIM):
        for a in range(DIM):
            for b in range(DIM):
                chris_contr[c] += hu[a][b] * chris[c][a][b]

    # Physical Christoffel symbols: computed directly as the Christoffel of the
    # physical spatial metric gamma_{ij} = h_{ij} / chi, from its definition
    #   gamma^{ij}    = chi h^{ij}
    #   d_l gamma_ij  = d_l h_ij / chi - h_ij d_l chi / chi^2
    #   Gamma^i_{jk}  = 1/2 gamma^{im} (d_j gamma_km + d_k gamma_jm - d_m gamma_jk)
    gam_UU = [[chi * hu[i][j] for j in range(DIM)] for i in range(DIM)]
    d_gam = zeros(DIM, DIM, DIM)  # d_gam[l][i][j] = d_l gamma_ij
    for l in range(DIM):
        for i in range(DIM):
            for j in range(DIM):
                d_gam[l][i][j] = (dh[l][i][j] / chi
                                  - h[i][j] * dchi[l] / chi ** 2)
    phys_chris = zeros(DIM, DIM, DIM)
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                s = 0.0
                for m in range(DIM):
                    s += gam_UU[i][m] * (d_gam[j][k][m] + d_gam[k][j][m]
                                         - d_gam[m][j][k])
                phys_chris[i][j][k] = 0.5 * s

    # P[a][c][e] = h_{ab} Gamma^b_{cd} h^{de}  (= Gamma_{ac}{}^e)
    P = zeros(DIM, DIM, DIM)
    for a in range(DIM):
        for c in range(DIM):
            for e in range(DIM):
                s = 0.0
                for b in range(DIM):
                    for d in range(DIM):
                        s += h[a][b] * chris[b][c][d] * hu[d][e]
                P[a][c][e] = s

    # Conformal Ricci (Alcubierre 2.8.17)
    ric_bar = zeros(DIM, DIM)
    for i in range(DIM):
        for j in range(DIM):
            # -1/2 h^{lm} d_l d_m h_{ij}
            term_a = 0.0
            for a in range(DIM):
                for b in range(DIM):
                    term_a += hu[a][b] * d2h[b][a][i][j]
            term_a *= -0.5
            # h_{k(i} d_{j)} Gamma^k  (symmetrised)
            term_b = 0.0
            for k in range(DIM):
                term_b += dgamma[i][k] * h[k][j] + dgamma[j][k] * h[k][i]
            term_b *= 0.5
            # 1/2 Gamma^m d_m h_{ij}
            term_c = 0.0
            for m in range(DIM):
                term_c += gamma[m] * dh[m][i][j]
            term_c *= 0.5
            # Gamma-Gamma terms
            r_ij = 0.0
            r_ji = 0.0
            for a in range(DIM):
                for b in range(DIM):
                    r_ij += chris[a][b][i] * P[j][a][b]
                    r_ji += chris[a][b][j] * P[i][a][b]
            term_d = r_ij + r_ji
            term_e = 0.0
            for a in range(DIM):
                for c in range(DIM):
                    term_e += chris[a][i][c] * P[a][j][c]
            ric_bar[i][j] = term_a + term_b + term_c + term_d + term_e

    # chi-dependent Ricci pieces
    box_chi = 0.0          # h^{ab} d_a d_b chi
    dchi_dchi = 0.0        # h^{ab} d_a chi d_b chi
    for a in range(DIM):
        for b in range(DIM):
            box_chi += hu[a][b] * d2chi[a][b]
            dchi_dchi += hu[a][b] * dchi[a] * dchi[b]
    scalar3 = (1.0 / (2 * chi) * box_chi
               - (dims - 1) / (4 * chi ** 2) * dchi_dchi)
    contr_dchi = sum(chris_contr[c] * dchi[c] for c in range(DIM))

    ric_chi = zeros(DIM, DIM)
    for i in range(DIM):
        for j in range(DIM):
            cov = d2chi[i][j]
            for m in range(DIM):
                cov -= dchi[m] * chris[m][i][j]
            ric_chi[i][j] = (
                (dims - 3) / (2 * chi) * cov
                - (dims - 3) / (4 * chi ** 2) * dchi[i] * dchi[j]
                + h[i][j] * scalar3
                - 1.0 / (2 * chi) * h[i][j] * contr_dchi)

    # Z-vector terms
    zz = [z_over_chi[k] * chi for k in range(DIM)]
    z_h = [sum(zz[k] * h[k][i] for k in range(DIM)) for i in range(DIM)]
    z_dchi = sum(zz[k] * dchi[k] for k in range(DIM))
    z_terms = zeros(DIM, DIM)
    for i in range(DIM):
        for j in range(DIM):
            z_terms[i][j] = chi ** (-2) * (
                z_h[i] * dchi[j] + z_h[j] * dchi[i] - h[i][j] * z_dchi)

    ricci_z = [[ric_bar[i][j] + ric_chi[i][j] + z_terms[i][j]
                for j in range(DIM)] for i in range(DIM)]
    ricci_z_scalar = chi * sum(hu[a][b] * ricci_z[a][b]
                               for a in range(DIM) for b in range(DIM))

    return dict(h_UU=hu, chris=chris, chris_contracted=chris_contr,
                phys_chris=phys_chris, ricciZ=ricci_z,
                ricciZ_scalar=ricci_z_scalar, z_terms=z_terms)


def d1_contracted_christoffel(h, hu, chris, dh, d2h):
    """d1cc[i][l] = d_l Gammabar^i, where Gammabar^i = h^{jk} Gamma^i_{jk} is the
    contracted conformal Christoffel.

    Derived by differentiating the definition (product + chain rule):
      d_l h^{jk}       = -h^{ja} h^{kb} d_l h_{ab}
      d_l Gamma^i_{jk} = 1/2 [ (d_l h^{im}) (d_j h_km + d_k h_jm - d_m h_jk)
                               + h^{im} (d_l d_j h_km + d_l d_k h_jm
                                         - d_l d_m h_jk) ]
      d_l Gammabar^i   = (d_l h^{jk}) Gamma^i_{jk} + h^{jk} (d_l Gamma^i_{jk})
    """
    # derivative of the inverse metric: duh[l][j][k] = d_l h^{jk}
    duh = zeros(DIM, DIM, DIM)
    for l in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                s = 0.0
                for a in range(DIM):
                    for b in range(DIM):
                        s += hu[j][a] * hu[k][b] * dh[l][a][b]
                duh[l][j][k] = -s

    # derivative of the conformal Christoffel: dchris[l][i][j][k] = d_l Gamma^i_jk
    dchris = zeros(DIM, DIM, DIM, DIM)
    for l in range(DIM):
        for i in range(DIM):
            for j in range(DIM):
                for k in range(DIM):
                    s = 0.0
                    for m in range(DIM):
                        s += duh[l][i][m] * (dh[j][k][m] + dh[k][j][m]
                                             - dh[m][j][k])
                        s += hu[i][m] * (d2h[l][j][k][m] + d2h[l][k][j][m]
                                         - d2h[l][m][j][k])
                    dchris[l][i][j][k] = 0.5 * s

    d1cc = zeros(DIM, DIM)  # d1cc[i][l] = d_l Gammabar^i
    for i in range(DIM):
        for l in range(DIM):
            s = 0.0
            for j in range(DIM):
                for k in range(DIM):
                    s += (duh[l][j][k] * chris[i][j][k]
                          + hu[j][k] * dchris[l][i][j][k])
            d1cc[i][l] = s
    return d1cc


def compute_ricci_general(v, e, dZ_coeff):
    """General-form Ricci Z, built independently on top of the already computed
    Ricci Z (e["ricciZ"]).

    Physically this is the conformal Ricci of the *metric* connection plus
    dZ_coeff times the symmetrised covariant derivative of the constraint vector
    Z_i = 1/2 (Gamma_i - Gammabar_i):

        R^gen_ij = R[Gammabar]_ij + dZ_coeff * D_(i Z_j)

    Here Gammabar^i = h^{jk} Gamma^i_jk is the metric-derived contracted
    connection and Gamma^i is the independent evolved connection variable.

    Building blocks (all independent of the code under test):
      * base = e["ricciZ"] - e["z_terms"] is the Ricci with the independent
        Gamma (compute_ricci_Z is linear in its Z vector, so subtracting its
        z_terms leaves the Gamma/chi contributions).
      * R[Gammabar] = base + 1/2 corr, where corr swaps the independent Gamma
        for the metric Gammabar in the connection terms of the Ricci.
      * 2 D_(i Z_j) = (base - R[Gammabar]) + z_terms(Z_constraint)/chi.
    """
    chi = v["chi"]
    h = v["h"]
    gamma = v["gamma"]
    dchi = v["dchi"]
    dh = v["dh"]
    dgamma = v["dgamma"]
    d2h = v["d2h"]
    hu = e["h_UU"]
    chris = e["chris"]
    gammabar = e["chris_contracted"]

    # Ricci with the independent Gamma: undo the random-Z terms in ricciZ
    base = [[e["ricciZ"][i][j] - e["z_terms"][i][j] for j in range(DIM)]
            for i in range(DIM)]

    # d_l Gammabar^i, derived from its definition
    d1cc = d1_contracted_christoffel(h, hu, chris, dh, d2h)

    # R[Gammabar]: swap Gamma -> Gammabar in the connection terms of the Ricci
    #   corr_ij = h_mi (d_j Gammabar^m - d_j Gamma^m)
    #           + h_mj (d_i Gammabar^m - d_i Gamma^m)
    #           + (Gammabar^m - Gamma^m) d_m h_ij
    r_metric = zeros(DIM, DIM)
    for i in range(DIM):
        for j in range(DIM):
            corr = 0.0
            for m in range(DIM):
                # d_j Gamma^m = dgamma[j][m]; d_m h_ij = dh[m][i][j]
                corr += (h[m][i] * (d1cc[m][j] - dgamma[j][m])
                         + h[m][j] * (d1cc[m][i] - dgamma[i][m])
                         + (gammabar[m] - gamma[m]) * dh[m][i][j])
            r_metric[i][j] = base[i][j] + 0.5 * corr

    # Constraint Z vector and the chi-part of its symmetrised covariant
    # derivative; assemble R[Gammabar] + dZ_coeff * D_(i Z_j).
    z_con = [0.5 * (gamma[i] - gammabar[i]) for i in range(DIM)]
    LL = zeros(DIM, DIM)
    for i in range(DIM):
        for j in range(DIM):
            zt = 0.0
            for k in range(DIM):
                zt += z_con[k] * (h[i][k] * dchi[j] + h[j][k] * dchi[i]
                                  - h[i][j] * dchi[k])
            two_dz = (base[i][j] - r_metric[i][j]) + zt / chi
            LL[i][j] = r_metric[i][j] + 0.5 * dZ_coeff * two_dz

    scalar = chi * sum(hu[i][j] * LL[i][j]
                       for i in range(DIM) for j in range(DIM))
    return LL, scalar


# --------------------------------------------------------------------------- #
# Header emission
# --------------------------------------------------------------------------- #
def num(x):
    """Full-precision, round-trippable C++ double literal."""
    return repr(float(x))


def emit_input_section(v):
    L = ["#ifdef CCZ4_GEOMETRY_INPUT_VALUES", ""]
    L.append("chi = {};".format(num(v["chi"])))
    L.append("")
    for i in range(DIM):
        for j in range(DIM):
            L.append("h({}, {}) = {};".format(i, j, num(v["h"][i][j])))
    L.append("")
    for i in range(DIM):
        L.append("Gamma({}) = {};".format(i, num(v["gamma"][i])))
    L.append("")
    for i in range(DIM):
        L.append("d1_chi({}) = {};".format(i, num(v["dchi"][i])))
    L.append("")
    L.append("// d1_h is stored symmetrically")
    for p, (i, j) in enumerate(SYM_PAIRS):
        for k in range(DIM):
            L.append("d1_h({}, {}) = {};".format(p, k, num(v["dh"][k][i][j])))
    L.append("")
    # d1_Gamma(a, b) = d_b Gamma^a ; internal dgamma[i][j] = d_i Gamma^j
    for a in range(DIM):
        for b in range(DIM):
            L.append("d1_Gamma({}, {}) = {};".format(
                a, b, num(v["dgamma"][b][a])))
    L.append("")

    for i in range(DIM):
        for j in range(DIM):
            L.append("d2_chi({}, {}) = {};".format(i, j, num(v["d2chi"][i][j])))
    L.append("")

    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    L.append("d2_h({}, {}, {}, {}) = {};".format(i, j, k, l, num(v["d2h"][k][l][i][j])))
    L.append("")

    for i in range(DIM):
        L.append("Z_over_chi({}) = {};".format(i, num(v["z_over_chi"][i])))
    L.append("")
    L.append("#endif // CCZ4_GEOMETRY_INPUT_VALUES")
    return L


def _row(vec):
    return "{" + ", ".join(num(x) for x in vec) + "}"


def _emit_rank2(L, decl, mat):
    L.append(decl + " = {")
    for i in range(DIM):
        comma = "," if i < DIM - 1 else ""
        L.append("    {}{}".format(_row(mat[i]), comma))
    L.append("};")


def _emit_rank3(L, decl, tens):
    L.append(decl + " = {")
    for i in range(DIM):
        for j in range(DIM):
            open_ = "    {" if j == 0 else "     "
            close = "}" if j == DIM - 1 else ""
            comma = "," if (j < DIM - 1 or i < DIM - 1) else ""
            L.append("{}{}{}{}".format(open_, _row(tens[i][j]), close, comma))
    L.append("};")


def emit_expected_section(e):
    L = ["#ifdef CCZ4_GEOMETRY_EXPECTED_VALUES", ""]

    _emit_rank2(L, "double h_UU_known[3][3]", e["h_UU"])
    _emit_rank3(L, "double chris_known[3][3][3]", e["chris"])

    L.append("double chris_contracted_known[3] = {};".format(
        _row(e["chris_contracted"])))

    _emit_rank3(L, "double chris_phys_known[3][3][3]", e["phys_chris"])

    _emit_rank2(L, "double ricciZ_known[3][3]", e["ricciZ"])

    L.append("double ricciZ_scalar_known = {};".format(num(e["ricciZ_scalar"])))
    L.append("")

    L.append("// General-form Ricci Z (compute_ricci_Z_general) with "
             "dZ_coeff = {:d}".format(e["dZ_coeff"]))
    _emit_rank2(L, "double ricciZ_general_known[3][3]", e["ricciZ_general"])
    L.append("double ricciZ_general_scalar_known = {};".format(
        num(e["ricciZ_general_scalar"])))
    L.append("")
    L.append("#endif // CCZ4_GEOMETRY_EXPECTED_VALUES")
    return L


HEADER = """\
/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Autogenerated by CCZ4GeometryGenerateExpectedValues.py
//
// This single header replaces CCZ4GeometryMathematicaValues.hpp and
// CCZ4GeometryMathematicaExpectedValues.hpp. It is deliberately NOT guarded
// against multiple inclusion: it is #included twice, once with each of the
// macros below defined, so it can be used at the two different points in the
// unit test.
//
//   #define CCZ4_GEOMETRY_INPUT_VALUES
//   #include "CCZ4GeometryExpectedValues.hpp"
//   #undef  CCZ4_GEOMETRY_INPUT_VALUES
//
//   #define CCZ4_GEOMETRY_EXPECTED_VALUES
//   #include "CCZ4GeometryExpectedValues.hpp"
//   #undef  CCZ4_GEOMETRY_EXPECTED_VALUES
//
// NOLINTBEGIN
"""

FOOTER = "// NOLINTEND\n"


def build_header(inputs, expected):
    lines = [HEADER]
    lines.extend(emit_input_section(inputs))
    lines.append("")
    lines.extend(emit_expected_section(expected))
    lines.append("")
    lines.append(FOOTER)
    return "\n".join(lines)


def sanity_check(inputs, expected):
    """Cheap internal consistency checks."""
    h, hu = inputs["h"], expected["h_UU"]
    assert abs(det3(h) - 1.0) < 1e-12, "conformal metric is not unit determinant"
    for i in range(DIM):
        for j in range(DIM):
            s = sum(h[i][k] * hu[k][j] for k in range(DIM))
            assert abs(s - (1.0 if i == j else 0.0)) < 1e-12, "h_UU not inverse"
    r = expected["ricciZ"]
    for i in range(DIM):
        for j in range(DIM):
            assert abs(r[i][j] - r[j][i]) < 1e-9, "ricciZ not symmetric"
    # Physical Christoffels: symmetric in the lower indices, and metric
    # compatible with gamma_{ij} = h_{ij}/chi, i.e.
    #   d_l gamma_ij = Gamma^m_{li} gamma_mj + Gamma^m_{lj} gamma_im.
    chi, dh, dchi = inputs["chi"], inputs["dh"], inputs["dchi"]
    pc = expected["phys_chris"]
    gam = [[h[i][j] / chi for j in range(DIM)] for i in range(DIM)]
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                assert abs(pc[i][j][k] - pc[i][k][j]) < 1e-9, \
                    "phys_chris not symmetric in lower indices"
    for l in range(DIM):
        for i in range(DIM):
            for j in range(DIM):
                d_gam = dh[l][i][j] / chi - h[i][j] * dchi[l] / chi ** 2
                rhs = sum(pc[m][l][i] * gam[m][j] + pc[m][l][j] * gam[i][m]
                          for m in range(DIM))
                assert abs(d_gam - rhs) < 1e-9, \
                    "phys_chris not metric compatible"
    if "ricciZ_general" in expected:
        rg = expected["ricciZ_general"]
        for i in range(DIM):
            for j in range(DIM):
                assert abs(rg[i][j] - rg[j][i]) < 1e-9, \
                    "ricciZ_general not symmetric"


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=None,
                        help="RNG seed for reproducible output (default: "
                             "nondeterministic, like the notebook)")
    parser.add_argument("--output", default=os.path.join(
        here, "CCZ4GeometryExpectedValues.hpp"),
        help="output header path")
    parser.add_argument("--dz-coeff", type=int, choices=(0, 1), default=1,
                        help="dZ_coeff for the general-form Ricci Z: a 0/1 bool "
                             "(0 -> pure metric Ricci, 1 -> full D_(i Z_j) "
                             "contribution); default 1")
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    inputs = generate_inputs()
    expected = compute_expected(inputs)
    ll_gen, scalar_gen = compute_ricci_general(inputs, expected, args.dz_coeff)
    expected["ricciZ_general"] = ll_gen
    expected["ricciZ_general_scalar"] = scalar_gen
    expected["dZ_coeff"] = args.dz_coeff
    sanity_check(inputs, expected)

    with open(args.output, "w") as f:
        f.write(build_header(inputs, expected))
    print("Wrote {}".format(args.output))


if __name__ == "__main__":
    main()
