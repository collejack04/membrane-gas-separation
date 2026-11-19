"""
Membrane gas separation — single-stage, binary, complete-mixing model.
Reference implementation for the GitHub Pages demo.

Fast gas A (e.g. CO2), slow gas B (e.g. N2 or CH4).
Feed: F, xF = mole fraction A in feed.
Permeate: P, y = mole fraction A in permeate.
Retentate: R, x = mole fraction A in retentate.
Stage cut: theta = P / F  (fraction of feed that permeates).

Core idea (the one line to remember):
    flux_i = permeance_i * driving force_i
    J_i = Q_i * (p_feed * x_i - p_perm * y_i)

Assumptions (complete-mixing / Weller-Steiner model, Seader et al.):
  1. Binary mixture, isothermal, steady state.
  2. Both sides well mixed, so retentate composition = x everywhere
     on the feed side, permeate composition = y everywhere.
  3. Permeance Q_i is constant (no plasticisation, no concentration
     polarisation, no pressure dependence).
  4. No pressure drop along the membrane.
  5. Selectivity alpha = Q_A / Q_B = P_A / P_B is constant.

Units:
  1 GPU = 1e-6 cm3(STP) / (cm2 s cmHg)
        = 3.35e-10 mol / (m2 s Pa)
  1 Barrer = 1e-10 cm3(STP) cm / (cm2 s cmHg)
  permeance [GPU] = permeability [Barrer] / thickness [micron]
  (because 1 Barrer / 1 micron = 1 GPU)

The separation itself depends only on alpha and pressure ratio
r = p_perm / p_feed, not on absolute permeance. Permeance (and
therefore permeability / thickness) sets the membrane AREA needed,
not the purity. That is why the demo has separate sliders:
selectivity / pressure move the purity-recovery curve, permeability
moves the area.
"""

from __future__ import annotations

import math

GPU_TO_SI = 3.35e-10  # mol / (m2 s Pa) per GPU
BAR_TO_PA = 1.0e5


def _weller_steiner_residual(y: float, x: float, alpha: float, r: float) -> float:
    """f(y) = 0 at the solution. y/(1-y) - alpha*(x - r*y)/((1-x) - r*(1-y))."""
    num = x - r * y
    den = (1.0 - x) - r * (1.0 - y)
    if den <= 1e-12 or num <= 0:
        return float("inf") if num <= 0 else -float("inf")
    return y / (1.0 - y) - alpha * num / den


def solve_stage(
    x_feed: float,
    alpha: float,
    pressure_ratio: float,
    theta: float,
) -> tuple[float, float]:
    """
    Solve the complete-mixing model for one stage cut.

    Returns (y, x): permeate purity, retentate composition.
    Uses the overall balance x = (xF - theta*y) / (1 - theta)
    substituted into the Weller-Steiner flux-ratio equation,
    solved by bisection on y in (xF, min(1, xF/theta)).
    """
    if not (0.0 < x_feed < 1.0):
        raise ValueError("x_feed must be in (0, 1)")
    if alpha <= 1.0:
        raise ValueError("alpha must be > 1 for A to be the fast gas")
    if not (0.0 <= pressure_ratio < 1.0):
        raise ValueError("pressure_ratio = p_perm / p_feed must be in [0, 1)")
    if not (0.0 < theta < 1.0):
        raise ValueError("theta must be in (0, 1)")

    # y cannot exceed xF/theta (that would need x < 0, i.e. >100% recovery)
    y_hi = min(0.999999, x_feed / theta * 0.999999)
    y_lo = x_feed * 1.000001 + 1e-9
    if y_hi <= y_lo:
        y_hi = min(0.999999, y_lo + 1e-6)

    def f(y: float) -> float:
        x = (x_feed - theta * y) / (1.0 - theta)
        if x <= 1e-12 or x >= 1.0 - 1e-12:
            return float("inf") if y > y_lo else -float("inf")
        return _weller_steiner_residual(y, x, alpha, pressure_ratio)

    # Bisection: f(y_lo) is normally negative, f(y_hi) positive.
    # Scan for a sign change to be safe at extreme theta / r.
    n_scan = 200
    prev_y = y_lo
    prev_f = f(prev_y)
    bracket = None
    for i in range(1, n_scan + 1):
        yy = y_lo + (y_hi - y_lo) * i / n_scan
        ff = f(yy)
        if prev_f == 0:
            bracket = (prev_y, prev_y)
            break
        if ff == 0:
            bracket = (yy, yy)
            break
        if (prev_f < 0 < ff) or (prev_f > 0 > ff):
            bracket = (prev_y, yy)
            break
        prev_y, prev_f = yy, ff
    if bracket is None:
        # No sign change found: fall back to the endpoint closest to zero.
        # This happens when theta is so high the feed is almost exhausted.
        return (y_hi, max(0.0, (x_feed - theta * y_hi) / (1.0 - theta)))

    lo, hi = bracket
    flo = f(lo)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if fm == 0 or (hi - lo) < 1e-10:
            lo = hi = mid
            break
        if (flo < 0 < fm) or (flo > 0 > fm):
            hi = mid
        else:
            lo, flo = mid, fm
    y = 0.5 * (lo + hi)
    x = (x_feed - theta * y) / (1.0 - theta)
    return y, max(0.0, min(1.0, x))


def evaluate(
    x_feed: float = 0.15,
    alpha: float = 30.0,
    permeability_barrer: float = 100.0,
    thickness_um: float = 0.1,
    feed_pressure_bar: float = 5.0,
    perm_pressure_bar: float = 1.0,
    theta: float = 0.5,
    feed_flow_mol_s: float = 10.0,
) -> dict:
    """Evaluate one operating point. Returns all demo outputs."""
    r = perm_pressure_bar / feed_pressure_bar
    y, x = solve_stage(x_feed, alpha, r, theta)

    permeance_gpu = permeability_barrer / thickness_um
    q_a_si = permeance_gpu * GPU_TO_SI
    q_b_si = q_a_si / alpha

    p_f = feed_pressure_bar * BAR_TO_PA
    p_p = perm_pressure_bar * BAR_TO_PA
    j_a = q_a_si * (p_f * x - p_p * y)
    j_b = q_b_si * (p_f * (1 - x) - p_p * (1 - y))
    j_tot = j_a + j_b

    perm_flow = theta * feed_flow_mol_s
    area_m2 = perm_flow / j_tot if j_tot > 0 else float("inf")
    recovery = theta * y / x_feed if x_feed > 0 else 0.0

    return {
        "x_feed": x_feed,
        "alpha": alpha,
        "permeability_barrer": permeability_barrer,
        "thickness_um": thickness_um,
        "permeance_gpu": permeance_gpu,
        "feed_pressure_bar": feed_pressure_bar,
        "perm_pressure_bar": perm_pressure_bar,
        "pressure_ratio": r,
        "theta": theta,
        "feed_flow_mol_s": feed_flow_mol_s,
        "y_permeate_purity": y,
        "x_retentate": x,
        "recovery_A": recovery,
        "permeate_flow_mol_s": perm_flow,
        "retentate_flow_mol_s": feed_flow_mol_s - perm_flow,
        "flux_A_mol_m2_s": j_a,
        "flux_B_mol_m2_s": j_b,
        "flux_total_mol_m2_s": j_tot,
        "area_m2": area_m2,
    }


def tradeoff_curve(
    x_feed: float = 0.15,
    alpha: float = 30.0,
    feed_pressure_bar: float = 5.0,
    perm_pressure_bar: float = 1.0,
    n: int = 40,
) -> list[dict]:
    """Purity vs recovery as theta is swept 0.02 -> 0.95."""
    r = perm_pressure_bar / feed_pressure_bar
    out = []
    for i in range(n):
        theta = 0.02 + (0.95 - 0.02) * i / (n - 1)
        y, x = solve_stage(x_feed, alpha, r, theta)
        out.append(
            {
                "theta": theta,
                "y": y,
                "x": x,
                "recovery": theta * y / x_feed,
            }
        )
    return out


if __name__ == "__main__":
    # Default demo case: 15% CO2 in N2, alpha=30, 100 Barrer / 0.1 um
    # (= 1000 GPU), 5 bar feed, 1 bar permeate, 50% stage cut,
    # 10 mol/s feed. Thin selective layer is why permeance is high:
    # permeance = permeability / thickness.
    res = evaluate()
    print("Default case: 15% CO2, alpha=30, 100 Barrer, 5 bar -> 1 bar, theta=0.50")
    for k, v in res.items():
        print(f"  {k:24s} = {v:.4g}" if isinstance(v, float) else f"  {k:24s} = {v}")
    print("\nTradeoff (theta, purity y, recovery):")
    for p in tradeoff_curve(n=9):
        print(f"  theta={p['theta']:.2f}  y={p['y']:.3f}  recovery={p['recovery']:.3f}")
