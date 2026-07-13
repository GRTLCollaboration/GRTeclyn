"""2D axisymmetric spinning Q-ball (Q-torus) solver.

Stationary rotating ansatz::

    Phi(t, rho, z, phi) = f(rho, z) * exp(i * (m_az * phi - omega * t))

Substituting into the flat-space Klein--Gordon equation with sextic potential
V = 1/2 mu^2 |Phi|^2 - 1/4 lam |Phi|^4 + 1/6 mu6 |Phi|^6 gives the 2D
elliptic PDE for the real modulus f(rho, z)::

    d^2f/drho^2 + (1/rho) df/drho + d^2f/dz^2
        - (m_az^2 / rho^2 + kappa^2) f + lam f^3 - mu6 f^5 = 0

with kappa^2 = mass^2 - omega^2 > 0 (localization), and BCs:
    f(0, z) = 0             (axis regularity for m_az >= 1)
    f(rho_max, z) = 0       (far-field decay)
    f(rho, z_max) = 0       (far-field decay)
    df/dz(rho, 0) = 0       (equatorial symmetry, ground state)

Solved by a damped Newton iteration with an explicit sparse Jacobian (direct
SuperLU solve), wrapped in an omega-continuation loop: the solve starts in the
thick-wall regime (omega close to mass, where the soliton amplitude is small and
the problem nearly linear so Newton converges from the twisted-spherical seed),
then omega is stepped down to the target while re-using the previous solution as
the seed. This is the standard existence path for rotating soliton eigenstates.

Output: QBallTorusProfile dataclass with the 2D grid and interpolant, plus
file I/O (``write_torus_profile`` / ``read_torus_profile``) for piping into
GRTresna's 2D torus loader.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from numpy.typing import NDArray
from scipy.interpolate import RectBivariateSpline
from scipy.sparse.linalg import splu

from .qball_couplings import QBallCouplings
from .qball_ode import QBallRadialProfile, solve_qball_radial_profile


@dataclass(frozen=True)
class QBallTorusProfile:
    """Tabulated 2D spinning Q-ball profile f(rho, z) on the (rho>=0, z>=0) half-plane."""

    mass: float
    lam: float
    mu: float
    omega: float
    m_az: int
    f_max: float
    rho: NDArray[np.float64]  # shape (n_rho,)
    z: NDArray[np.float64]  # shape (n_z,)
    f: NDArray[np.float64]  # shape (n_rho, n_z)
    spline: RectBivariateSpline

    @property
    def noether_charge(self) -> float:
        """Integrated Noether charge Q = 2 omega int |Phi|^2 d^3x over full space.

        With |Phi| = f(rho,z) and the z>=0 half-plane stored, the full-space
        integral picks up a factor 2 (z-reflection) and 2 pi (azimuthal)::

            Q = 2 omega * (2 pi) * 2 * int_0^inf int_0^inf f^2 rho drho dz
        """
        drho = float(self.rho[1] - self.rho[0]) if len(self.rho) > 1 else 1.0
        dz = float(self.z[1] - self.z[0]) if len(self.z) > 1 else 1.0
        integrand = self.f**2 * self.rho[:, None]
        return 8.0 * math.pi * self.omega * float(np.sum(integrand)) * drho * dz

    @property
    def energy(self) -> float:
        """Field energy E = int T_00 d^3x with the same half-plane convention."""
        drho = float(self.rho[1] - self.rho[0]) if len(self.rho) > 1 else 1.0
        dz = float(self.z[1] - self.z[0]) if len(self.z) > 1 else 1.0
        f = self.f
        rho = self.rho[:, None]

        df_drho = np.zeros_like(f)
        df_dz = np.zeros_like(f)
        df_drho[1:-1, :] = (f[2:, :] - f[:-2, :]) / (2.0 * drho)
        df_dz[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * dz)

        kinetic = self.omega**2 * f**2
        gradient = (
            df_drho**2
            + df_dz**2
            + (self.m_az**2 / np.maximum(rho**2, 1e-30)) * f**2
        )
        potential = (
            0.5 * self.mass**2 * f**2
            - 0.25 * self.lam * f**4
            + (1.0 / 6.0) * self.mu * f**6
        )
        integrand = (kinetic + gradient + potential) * rho
        return 8.0 * math.pi * float(np.sum(integrand)) * drho * dz

    def eval_f(self, rho_q: float | NDArray, z_q: float | NDArray) -> NDArray[np.float64]:
        """Evaluate f at arbitrary (rho, z) via the bivariate spline (|z| symmetry)."""
        rho_arr = np.atleast_1d(np.asarray(rho_q, dtype=np.float64))
        z_arr = np.abs(np.atleast_1d(np.asarray(z_q, dtype=np.float64)))
        rho_arr = np.clip(rho_arr, float(self.rho[0]), float(self.rho[-1]))
        z_arr = np.clip(z_arr, float(self.z[0]), float(self.z[-1]))
        return np.asarray(self.spline(rho_arr, z_arr, grid=False), dtype=np.float64)


def _build_seed(
    radial: QBallRadialProfile,
    rho_int: NDArray[np.float64],
    z_int: NDArray[np.float64],
    m_az: int,
) -> NDArray[np.float64]:
    """Seed from the 1D spherical profile: f0 = phi_1D(r) * (rho/r)^m_az."""
    RHO, Z = np.meshgrid(rho_int, z_int, indexing="ij")
    R = np.sqrt(RHO**2 + Z**2)
    R_safe = np.maximum(R, 1.0e-12)
    sin_theta = RHO / R_safe
    phi_1d = radial.eval_phi0(R.ravel()).reshape(R.shape)
    seed = phi_1d * sin_theta**m_az
    return np.maximum(seed, 0.0)


def _build_linear_operator(
    rho_int: NDArray[np.float64],
    n_rho: int,
    n_z: int,
    drho: float,
    dz: float,
    m_az: int,
) -> sp.csc_matrix:
    """Assemble the omega-independent sparse operator L0 = Laplacian - centrifugal.

    (The -kappa^2 term is added separately so omega can vary in the bordered solve.)

    Unknowns are the interior points f[i, j], i in [0, n_rho), j in [0, n_z),
    flattened as k = i * n_z + j.  Grid layout:
      rho index i -> rho_int[i]  (rho=0 axis and rho_max are Dirichlet ghosts)
      z index j=0 is the equator (Neumann df/dz=0); j=n_z is z_max (Dirichlet).
    """
    drho2 = drho * drho
    dz2 = dz * dz
    n = n_rho * n_z

    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []

    def add(k: int, kk: int, v: float) -> None:
        rows.append(k)
        cols.append(kk)
        vals.append(v)

    for i in range(n_rho):
        rho_i = rho_int[i]
        c_rho2 = 1.0 / drho2
        c_rho1 = 1.0 / (2.0 * rho_i * drho)
        coeff_ip = c_rho2 + c_rho1  # coefficient of f[i+1, j]
        coeff_im = c_rho2 - c_rho1  # coefficient of f[i-1, j]
        centrifugal = (m_az * m_az) / (rho_i * rho_i)
        for j in range(n_z):
            k = i * n_z + j
            diag = -2.0 / drho2 - 2.0 / dz2 - centrifugal

            # rho neighbours (Dirichlet 0 at i=-1 and i=n_rho -> just skip)
            if i + 1 < n_rho:
                add(k, (i + 1) * n_z + j, coeff_ip)
            if i - 1 >= 0:
                add(k, (i - 1) * n_z + j, coeff_im)

            # z neighbours
            if j + 1 < n_z:
                add(k, i * n_z + (j + 1), 1.0 / dz2)
            # j+1 == n_z is Dirichlet 0 -> skip
            if j - 1 >= 0:
                add(k, i * n_z + (j - 1), 1.0 / dz2)
            else:
                # j == 0: Neumann df/dz=0 -> ghost f[i,-1] = f[i,0], folds into diagonal
                diag += 1.0 / dz2

            add(k, k, diag)

    L0 = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
    return L0


def _linop_diag_kappa(
    L0: sp.csc_matrix, kappa_sq: float, n: int
) -> sp.csc_matrix:
    """Return L0 with the -kappa^2 term on the diagonal (L0 excludes kappa^2)."""
    return L0 - sp.eye(n, format="csc") * kappa_sq


def _newton_solve_amplitude(
    L0: sp.csc_matrix,
    f0: NDArray[np.float64],
    omega0: float,
    *,
    mass: float,
    lam: float,
    mu: float,
    amp: float,
    k_ref: int,
    tol: float,
    maxiter: int,
    verbose: bool,
) -> tuple[NDArray[np.float64], float, float, bool]:
    """Bordered Newton solving (f, omega) at fixed peak amplitude f[k_ref]=amp.

    L0 is the omega-independent operator (Laplacian - centrifugal). The full
    linear operator is L(omega) = L0 - kappa^2 I with kappa^2 = mass^2 - omega^2.

    Residual system (size n+1):
        R_i    = [L(omega) f + lam f^3 - mu f^5]_i          (i = 0..n-1)
        R_ampl = f[k_ref] - amp

    Jacobian (bordered):
        [ J           | dR/domega ]   with J = L(omega) + diag(3 lam f^2 - 5 mu f^4)
        [ e_kref^T    | 0         ]        dR/domega = 2 omega f   (since d(-kappa^2 f)/domega)

    Returns (f, omega, resnorm, converged).
    """
    f = f0.copy()
    omega = float(omega0)
    n = f.size
    resnorm = np.inf
    converged = False
    e_kref = np.zeros(n)
    e_kref[k_ref] = 1.0

    for it in range(maxiter):
        kappa_sq = mass * mass - omega * omega
        L = _linop_diag_kappa(L0, kappa_sq, n)
        res_pde = L @ f + lam * f**3 - mu * f**5
        res_amp = f[k_ref] - amp
        resnorm = max(float(np.max(np.abs(res_pde))), abs(res_amp))
        if verbose:
            print(f"    Newton {it}: ||R||={resnorm:.3e} omega={omega:.5f} f_ref={f[k_ref]:.5f}")
        if resnorm < tol:
            converged = True
            break

        diag_nl = 3.0 * lam * f**2 - 5.0 * mu * f**4
        J = L + sp.diags(diag_nl, 0, shape=(n, n), format="csc")
        b = (2.0 * omega * f).reshape(n, 1)  # dR_pde/domega
        # Bordered sparse system: [[J, b],[e_kref^T, 0]]
        J_aug = sp.bmat(
            [[J, sp.csc_matrix(b)], [sp.csc_matrix(e_kref.reshape(1, n)), None]],
            format="csc",
        )
        rhs = np.concatenate([-res_pde, [-res_amp]])
        try:
            delta = splu(J_aug).solve(rhs)
        except Exception:
            delta = splu(J_aug + sp.eye(n + 1, format="csc") * 1e-8).solve(rhs)

        df = delta[:n]
        domega = float(delta[n])

        # Damped line search on the augmented residual norm
        step = 1.0
        accepted = False
        for _ in range(30):
            f_try = f + step * df
            omega_try = omega + step * domega
            if omega_try <= 0.0 or omega_try >= mass:
                step *= 0.5
                continue
            ks = mass * mass - omega_try * omega_try
            L_try = _linop_diag_kappa(L0, ks, n)
            rp = L_try @ f_try + lam * f_try**3 - mu * f_try**5
            ra = f_try[k_ref] - amp
            nrm = max(float(np.max(np.abs(rp))), abs(ra))
            if nrm < resnorm:
                f, omega = f_try, omega_try
                accepted = True
                break
            step *= 0.5
        if not accepted:
            f = f + step * df
            omega = min(max(omega + step * domega, 1e-6), mass - 1e-6)
    return f, omega, resnorm, converged


def solve_qball_torus_profile(
    couplings: QBallCouplings,
    *,
    m_az: int = 1,
    n_rho: int = 160,
    n_z: int = 160,
    rho_max: float | None = None,
    z_max: float | None = None,
    tol: float = 1.0e-7,
    maxiter: int = 60,
    omega_tol: float = 2.0e-3,
    max_amp_steps: int = 40,
    verbose: bool = False,
) -> QBallTorusProfile:
    """Solve for the stationary rotating Q-torus eigenstate f(rho, z).

    Method: **amplitude continuation**. For a fixed peak amplitude ``amp`` the
    bordered Newton system solves jointly for (f, omega); this cannot collapse to
    the trivial vacuum. Starting from a small amplitude (near the vacuum
    bifurcation, omega ~ mass) the amplitude is increased geometrically, lowering
    omega, until the target ``couplings.omega`` is bracketed; a final bisection on
    ``amp`` lands on the target omega within ``omega_tol``.
    """
    if m_az < 1:
        raise ValueError("m_az must be >= 1 for a spinning Q-torus")
    couplings.validate()
    if couplings.lam <= 0.0 or couplings.mu <= 0.0:
        raise ValueError("Q-torus solver requires lam > 0 and mu > 0")

    mass = couplings.mass
    lam = couplings.lam
    mu = couplings.mu
    omega_target = couplings.omega
    if omega_target >= mass:
        raise ValueError("omega must be < mass for localization")

    width = couplings.bound_state_width
    if rho_max is None:
        rho_max = max(30.0, 12.0 * width)
    if z_max is None:
        z_max = max(30.0, 12.0 * width)

    # Interior grids: rho in (0, rho_max) exclusive; z in [0, z_max) with z=0 equator
    rho_full = np.linspace(0.0, rho_max, n_rho + 2)
    rho_int = rho_full[1:-1].copy()
    drho = float(rho_full[1] - rho_full[0])
    z_full = np.linspace(0.0, z_max, n_z + 1)
    z_int = z_full[:-1].copy()
    dz = float(z_full[1] - z_full[0])

    L0 = _build_linear_operator(rho_int, n_rho, n_z, drho, dz, m_az)
    n = n_rho * n_z

    # Analytic toroidal seed peaked off-axis at rho0 ~ few * width (m=1 vanishes
    # on the axis, hence the (rho/rho0) prefactor).  The amplitude pin is placed
    # at this off-axis peak and held fixed through the continuation.
    core = couplings.core_amplitude if couplings.core_amplitude > 0 else 0.1
    rho0 = max(2.0 * width, 3.0)
    sig = max(1.5 * width, 2.0)
    RHO, Z = np.meshgrid(rho_int, z_int, indexing="ij")
    seed = (RHO / rho0) * np.exp(-((RHO - rho0) ** 2 + Z**2) / (2.0 * sig * sig))
    seed = np.maximum(seed, 0.0)
    seed /= max(seed.max(), 1e-30)
    k_ref = int(np.argmax(seed.ravel()))
    i_peak, j_peak = divmod(k_ref, n_z)

    if verbose:
        print(
            f"Grid: {n_rho}x{n_z} interior, drho={drho:.4f}, dz={dz:.4f}, "
            f"rho_max={rho_max:.1f}, z_max={z_max:.1f}"
        )
        print(f"seed torus rho0={rho0:.2f} sig={sig:.2f}; core_amp={core:.5f}")
        print(f"amplitude pin at rho={rho_int[i_peak]:.3f}, z={z_int[j_peak]:.3f}")
        print(f"target omega={omega_target:.4f}, mass={mass}")

    def boundary_contact(fvec: NDArray[np.float64]) -> float:
        """Max |f| on the outer rho / z edges relative to peak (soliton escapes box)."""
        g = fvec.reshape(n_rho, n_z)
        edge = max(np.max(np.abs(g[-1, :])), np.max(np.abs(g[:, -1])))
        return float(edge / max(np.max(np.abs(g)), 1e-30))

    # Amplitude continuation: start at a moderate fraction of the core amplitude
    # and grow.  omega decreases monotonically with amplitude; stop when the
    # target omega is bracketed, then bisect on amplitude.
    amp = 0.4 * core
    f = seed.ravel() * amp
    omega = 0.5 * mass

    hist: list[tuple[float, float, NDArray[np.float64]]] = []
    bracket: tuple | None = None
    for _ in range(max_amp_steps):
        f, omega, resnorm, conv = _newton_solve_amplitude(
            L0, f, omega, mass=mass, lam=lam, mu=mu, amp=amp,
            k_ref=k_ref, tol=tol, maxiter=maxiter, verbose=False,
        )
        f = np.maximum(f, 0.0)
        bc = boundary_contact(f)
        if verbose:
            print(f"  amp={amp:.5f} -> omega={omega:.5f} "
                  f"||R||={resnorm:.2e} conv={conv} f_max={f.max():.5f} edge={bc:.1e}")
        if not conv:
            raise RuntimeError(
                f"Q-torus bordered Newton failed at amp={amp:.5f}: ||R||={resnorm:.3e}"
            )
        if bc > 1.0e-2:
            raise RuntimeError(
                f"Q-torus escaping the box (edge/peak={bc:.2e}) at omega={omega:.4f}; "
                f"omega={omega_target:.4f} is too deep thin-wall — raise omega_target "
                f"(more compact) or increase rho_max/z_max."
            )
        hist.append((amp, omega, f.copy()))
        if abs(omega - omega_target) < omega_tol:
            break
        if omega < omega_target and len(hist) >= 2:
            bracket = (hist[-2], hist[-1])
            break
        amp = min(amp * 1.15, 0.999 * core)
        if amp >= 0.999 * core and hist and hist[-1][0] >= 0.999 * core:
            raise RuntimeError(
                f"Q-torus: hit core-amplitude ceiling (omega={omega:.4f}) before reaching "
                f"target {omega_target:.4f}; target is below omega_min-like floor for "
                f"these couplings — raise omega_target."
            )

    if bracket is not None and abs(omega - omega_target) >= omega_tol:
        (amp_lo, om_lo, f_lo), (amp_hi, om_hi, f_hi) = bracket
        for _ in range(40):
            amp_mid = 0.5 * (amp_lo + amp_hi)
            f_seed = f_hi if abs(om_hi - omega_target) < abs(om_lo - omega_target) else f_lo
            f, omega, resnorm, conv = _newton_solve_amplitude(
                L0, f_seed.copy(), 0.5 * (om_lo + om_hi), mass=mass, lam=lam, mu=mu,
                amp=amp_mid, k_ref=k_ref, tol=tol, maxiter=maxiter, verbose=False,
            )
            f = np.maximum(f, 0.0)
            if verbose:
                print(f"  bisect amp={amp_mid:.5f} -> omega={omega:.5f}")
            if abs(omega - omega_target) < omega_tol:
                break
            if omega > omega_target:
                amp_lo, om_lo, f_lo = amp_mid, omega, f.copy()
            else:
                amp_hi, om_hi, f_hi = amp_mid, omega, f.copy()

    if abs(omega - omega_target) >= 5.0 * omega_tol:
        raise RuntimeError(
            f"Q-torus: could not reach target omega={omega_target:.4f} "
            f"(got {omega:.4f}); adjust omega_target / grid."
        )

    f_2d = f.reshape(n_rho, n_z)
    f_full = np.zeros((n_rho + 2, n_z + 1), dtype=np.float64)
    f_full[1:-1, :-1] = f_2d

    spline = RectBivariateSpline(rho_full, z_full, f_full, kx=3, ky=3)
    f_max = float(np.max(f_full))

    profile = QBallTorusProfile(
        mass=mass,
        lam=lam,
        mu=mu,
        omega=float(omega),
        m_az=m_az,
        f_max=f_max,
        rho=rho_full,
        z=z_full,
        f=f_full,
        spline=spline,
    )
    if verbose:
        print(f"Converged. omega={omega:.5f} f_max={f_max:.6f}  "
              f"Q={profile.noether_charge:.4f}  E={profile.energy:.4f}")
    return profile


# ---------------------------------------------------------------------------
# File I/O — interop with GRTresna torus loader
# ---------------------------------------------------------------------------

_TORUS_HEADER = "# qball_torus_2d"


def write_torus_profile(profile: QBallTorusProfile, path: str | Path) -> None:
    """Write a 2D torus profile in the format expected by GRTresna.

    Layout: metadata comment header, then ``n_rho*n_z`` rows of ``rho z f`` with
    z varying fastest (row-major over (i_rho, j_z)).
    """
    path = Path(path)
    n_rho = len(profile.rho)
    n_z = len(profile.z)
    with open(path, "w") as fh:
        fh.write(f"{_TORUS_HEADER}\n")
        fh.write(
            f"# mass={profile.mass} lam={profile.lam} mu6={profile.mu} "
            f"omega={profile.omega} m_az={profile.m_az}\n"
        )
        fh.write(
            f"# n_rho={n_rho} n_z={n_z} "
            f"rho_max={float(profile.rho[-1]):.8g} z_max={float(profile.z[-1]):.8g}\n"
        )
        fh.write(
            f"# f_max={profile.f_max:.8g} Q={profile.noether_charge:.8g} "
            f"E={profile.energy:.8g}\n"
        )
        for i in range(n_rho):
            for j in range(n_z):
                fh.write(
                    f"{profile.rho[i]:.8g} {profile.z[j]:.8g} {profile.f[i, j]:.10g}\n"
                )


def read_torus_profile(path: str | Path) -> QBallTorusProfile:
    """Read a 2D torus profile written by :func:`write_torus_profile`."""
    path = Path(path)
    meta: dict[str, str] = {}
    data_lines: list[str] = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("#"):
                for token in s.lstrip("# ").split():
                    if "=" in token:
                        key, val = token.split("=", 1)
                        meta[key] = val
            elif s:
                data_lines.append(s)

    n_rho = int(meta["n_rho"])
    n_z = int(meta["n_z"])
    raw = np.loadtxt(data_lines, dtype=np.float64)
    if raw.shape != (n_rho * n_z, 3):
        raise ValueError(f"Expected {n_rho*n_z} rows x3, got {raw.shape}")

    rho_vals = raw[::n_z, 0].copy()
    z_vals = raw[:n_z, 1].copy()
    f_vals = raw[:, 2].reshape(n_rho, n_z).copy()
    spline = RectBivariateSpline(rho_vals, z_vals, f_vals, kx=3, ky=3)

    return QBallTorusProfile(
        mass=float(meta["mass"]),
        lam=float(meta["lam"]),
        mu=float(meta["mu6"]),
        omega=float(meta["omega"]),
        m_az=int(meta["m_az"]),
        f_max=float(np.max(f_vals)),
        rho=rho_vals,
        z=z_vals,
        f=f_vals,
        spline=spline,
    )
