#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mass–spring parameter estimation (super simple)
Methods:
  A) Derivative regression on: ydd = -a*yd - b*y  ->  omega0 = sqrt(b), zeta = a/(2*omega0)
  B) Nonlinear fit by integrating the ODE: m*ydd + b*yd + k*y = 0 (fit k, b)
Usage:
  - As a library: import mass_spring_fit and call fit_derivative_regression(...) or fit_ode(...)
  - As a script:  python3 mass_spring_fit.py data.csv --mass 0.200
CSV format:
  Expect columns named 't' and 'x' (in seconds and meters). Extra columns are ignored.
  If you don't have a CSV yet, run with --demo to generate synthetic data and fit it.
"""
import argparse
from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

@dataclass
class FitResult:
    omega0: float
    zeta: float
    k: float
    b_visc: float
    y0: float
    v0: float

def _odd(n: int) -> int:
    """Return nearest odd >= 5 (for Savitzky–Golay)."""
    n = max(5, n)
    return n if n % 2 == 1 else n + 1

def _auto_window_length(N: int) -> int:
    # Aim ~ 1/10 of samples, odd, at least 5
    wl = max(5, (N // 10) | 1)  # make it odd
    # keep window < N and reasonable upper bound
    wl = min(wl, N - (1 - (N % 2)))  # ensure < N and odd feasibility
    if wl < 5:
        wl = _odd(5)
    if wl >= N:
        wl = _odd(N - 1) if N > 5 else 5
    return int(wl)

def center_series(y: np.ndarray) -> np.ndarray:
    return y - np.median(y)

def fit_derivative_regression(t: np.ndarray, x: np.ndarray, poly: int = 3, window_length: int | None = None) -> FitResult:
    """
    Method A: Regress ydd = -a*yd - b*y using Savitzky–Golay derivatives.
    Returns omega0, zeta, and (if mass is provided later) k, b_visc can be derived.
    """
    t = np.asarray(t).astype(float)
    x = np.asarray(x).astype(float)
    y = center_series(x)
    # uniform dt assumption
    dt = np.median(np.diff(t))
    if dt <= 0:
        raise ValueError("Non-increasing time vector")
    N = len(t)
    if window_length is None:
        window_length = _auto_window_length(N)
    window_length = int(window_length)
    if window_length >= N:
        window_length = _odd(N - 1) if N > 5 else 5
    # Smooth and derivatives
    y_s = savgol_filter(y, window_length=window_length, polyorder=min(poly, window_length-1))
    yd = savgol_filter(y, window_length=window_length, polyorder=min(poly, window_length-1), deriv=1, delta=dt)
    ydd = savgol_filter(y, window_length=window_length, polyorder=min(poly, window_length-1), deriv=2, delta=dt)
    # Build regression: ydd = -a*yd - b*y  ->  X*[a,b]^T = -ydd
    X = np.column_stack([yd, y_s])
    rhs = -ydd
    # Solve least squares
    theta, *_ = np.linalg.lstsq(X, rhs, rcond=None)
    a, b = theta
    # Map to physical
    omega0 = np.sqrt(abs(b))
    if omega0 == 0:
        zeta = 0.0
    else:
        zeta = a / (2.0 * omega0)
    # Initial conditions (rough)
    y0 = y[0]
    v0 = (y[1] - y[0]) / dt if N > 1 else 0.0
    # k and b_visc undefined without mass; fill with NaN here
    return FitResult(omega0=float(omega0), zeta=float(zeta), k=np.nan, b_visc=np.nan, y0=float(y0), v0=float(v0))

def _simulate_ode(t: np.ndarray, y0: float, v0: float, m: float, k: float, b: float) -> np.ndarray:
    """Integrate m*y'' + b*y' + k*y = 0 over given t, return y(t)."""
    def f(_t, z):
        y, v = z
        return [v, -(b/m)*v - (k/m)*y]
    sol = solve_ivp(f, (t[0], t[-1]), [y0, v0], t_eval=t, method="RK45", rtol=1e-6, atol=1e-9)
    if not sol.success:
        raise RuntimeError(f"ODE solver failed: {sol.message}")
    return sol.y[0]

def fit_ode(t: np.ndarray, x: np.ndarray, m: float, guess: tuple[float, float] | None = None) -> FitResult:
    """
    Method B: Fit k,b by integrating the ODE and minimizing squared residuals.
    m: effective mass (kg). Provide m ≈ mass + m_spring/3.
    guess: optional (k0, b0). If None, derive from method A.
    """
    t = np.asarray(t).astype(float)
    x = np.asarray(x).astype(float)
    y = center_series(x)
    dt = np.median(np.diff(t))
    if dt <= 0:
        raise ValueError("Non-increasing time vector")
    # Initial conditions
    y0 = y[0]
    v0 = (y[1] - y[0]) / dt if len(t) > 1 else 0.0
    # Initial guess from method A
    if guess is None:
        fa = fit_derivative_regression(t, y)
        omega0 = max(1e-6, fa.omega0)
        zeta = max(0.0, fa.zeta)
        k0 = m * omega0**2
        b0 = 2 * m * zeta * omega0
    else:
        k0, b0 = guess
    # Bounds: keep positive k,b
    lb = np.array([1e-9, 0.0])
    ub = np.array([np.inf, np.inf])
    def residuals(theta):
        k, b = theta
        yhat = _simulate_ode(t, y0, v0, m, k, b)
        return yhat - y
    res = least_squares(residuals, x0=np.array([k0, b0]), bounds=(lb, ub), method="trf", loss="huber", f_scale=1.0, xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=200)
    k_fit, b_fit = res.x
    # Map to omega0, zeta
    omega0 = np.sqrt(k_fit / m)
    zeta = b_fit / (2.0 * m * omega0) if omega0 > 0 else 0.0
    return FitResult(omega0=float(omega0), zeta=float(zeta), k=float(k_fit), b_visc=float(b_fit), y0=float(y0), v0=float(v0))

def _load_csv(path: str):
    df = pd.read_csv(path)
    # try to find columns t and x (case-insensitive)
    cols = {c.lower(): c for c in df.columns}
    if 't' in cols and 'x' in cols:
        t = df[cols['t']].to_numpy(dtype=float)
        x = df[cols['x']].to_numpy(dtype=float)
    else:
        # fallback to first two numeric columns
        num_cols = [c for c in df.columns if np.issubdtype(df[c].dtype, np.number)]
        if len(num_cols) < 2:
            raise ValueError("CSV must contain columns 't' and 'x' or at least two numeric columns")
        t = df[num_cols[0]].to_numpy(dtype=float)
        x = df[num_cols[1]].to_numpy(dtype=float)
    return t, x

def _demo_data(T=5.0, fs=200, m=0.200, k=10.0, zeta=0.03, A=0.05, phi=0.2, noise=0.001, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, T, 1.0/fs)
    omega0 = np.sqrt(k/m)
    gamma = zeta * 2 * omega0
    omega_d = np.sqrt(max(omega0**2 - (gamma/2)**2, 1e-12))
    y = A * np.exp(-gamma * t/2) * np.cos(omega_d * t + phi)
    x = y + rng.normal(0.0, noise, size=t.size) + 0.0  # add tiny noise
    return t, x

def main():
    ap = argparse.ArgumentParser(description="Fit mass–spring parameters from time series x(t).")
    ap.add_argument("csv", nargs="?", help="CSV with columns t,x (seconds, meters).")
    ap.add_argument("--mass", "-m", type=float, default=0.200, help="Effective mass in kg (default: 0.200)")
    ap.add_argument("--demo", action="store_true", help="Run a synthetic demo instead of loading CSV")
    ap.add_argument("--window", type=int, default=None, help="Savitzky–Golay window length (odd). If omitted, auto.")
    args = ap.parse_args()

    if args.demo or not args.csv:
        print("[demo] Generating synthetic dataset...")
        t, x = _demo_data()
    else:
        t, x = _load_csv(args.csv)

    # Method A
    Ares = fit_derivative_regression(t, x, window_length=args.window)
    kA = args.mass * (Ares.omega0 ** 2)
    bA = 2 * args.mass * Ares.zeta * Ares.omega0
    print("\n=== Method A: Derivative regression ===")
    print(f"omega0  ≈ {Ares.omega0:.6f} rad/s")
    print(f"zeta    ≈ {Ares.zeta:.6f}  (-)")
    print(f"k (A)   ≈ {kA:.6f} N/m")
    print(f"b (A)   ≈ {bA:.6f} N·s/m")

    # Method B
    Bres = fit_ode(t, x, m=args.mass, guess=(kA, bA))
    print("\n=== Method B: ODE fit (nonlinear) ===")
    print(f"omega0  ≈ {Bres.omega0:.6f} rad/s")
    print(f"zeta    ≈ {Bres.zeta:.6f}  (-)")
    print(f"k (B)   ≈ {Bres.k:.6f} N/m")
    print(f"b (B)   ≈ {Bres.b_visc:.6f} N·s/m")

if __name__ == "__main__":
    main()
