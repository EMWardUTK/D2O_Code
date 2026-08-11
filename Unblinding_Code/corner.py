import argparse
import os
import re
import sys
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Parameter order for the grid (top-left to bottom-right on the diagonal).
# Only those with a found 1D file are kept, in this order.
DEFAULT_ORDER = ["vD", "vO", "vFe", "vPb", "BRN", "SSB"]

# Pretty axis labels
LABELS = {
    "vD":  r"$N_{\nu D}$",
    "vO":  r"$N_{\nu O}$",
    "vFe": r"$N_{\nu Fe}$",
    "vPb": r"$N_{\nu Pb}$",
    "BRN": r"$N_{\mathrm{BRN}}$",
    "SSB": r"$N_{\mathrm{SSB}}$",
}

# Delta(2lnL) confidence levels
LEVELS_1D = [1.00, 3.84]            # 68.3%, 95% for 1 dof
LEVELS_2D = [2.30, 6.18]           # 68.3%, 95.4% for 2 dof

def interval_1d(x, y, thr):
    """Find threshold crossings of (y - y_min) via linear interpolation."""
    z = y - np.nanmin(y)
    crossings = []
    for i in range(len(x) - 1):
        a, b = z[i] - thr, z[i + 1] - thr
        if np.isnan(a) or np.isnan(b):
            continue
        if a == 0:
            crossings.append(x[i])
        if a * b < 0:
            t = a / (a - b)
            crossings.append(x[i] + t * (x[i + 1] - x[i]))
    return crossings

def load_1d(path):
    """Return (value, twoDeltaNLL) arrays sorted by value."""
    d = np.loadtxt(path, comments="#")
    d = d[np.argsort(d[:, 0])]
    return d[:, 0], d[:, 1]


def load_2d(path):
    """Return (X, Y, Z) grids suitable for plt.contour.

    The C++ writer loops x outer, y inner, so the file is grouped by x.
    We reshape onto a regular grid by reading unique sorted axis values.
    """
    d = np.loadtxt(path, comments="#")
    xs = np.unique(d[:, 0])
    ys = np.unique(d[:, 1])
    nx, ny = len(xs), len(ys)

    Z = np.full((ny, nx), np.nan)
    # map each value to its index for robustness (handles any ordering)
    xi = {v: i for i, v in enumerate(xs)}
    yi = {v: i for i, v in enumerate(ys)}
    for x, y, z in d:
        Z[yi[y], xi[x]] = z

    X, Y = np.meshgrid(xs, ys)
    return X, Y, Z


def discover(indir, order):
    """Find which parameters have 1D files; keep requested order."""
    present = []
    for name in order:
        if os.path.exists(os.path.join(indir, f"pll_1D_{name}.txt")):
            present.append(name)
    # also catch any 1D files not in the default order
    for fn in sorted(os.listdir(indir)):
        m = re.match(r"pll_1D_(.+)\.txt$", fn)
        if m and m.group(1) not in present:
            present.append(m.group(1))
    return present


def find_2d(indir, a, b):
    """Return (path, swapped) for the 2D file of pair (a, b), or (None, False).

    swapped=True means the file stores (b, a) and axes must be transposed.
    """
    p1 = os.path.join(indir, f"pll_2D_{a}_{b}.txt")
    p2 = os.path.join(indir, f"pll_2D_{b}_{a}.txt")
    if os.path.exists(p1):
        return p1, False
    if os.path.exists(p2):
        return p2, True
    return None, False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="pll_output", help="directory with pll_*.txt")
    ap.add_argument("--out", default="pll_corner.png", help="output image path")
    ap.add_argument("--order", default=",".join(DEFAULT_ORDER),
                    help="comma-separated parameter order")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()

    if not os.path.isdir(args.indir):
        sys.exit(f"Input directory not found: {args.indir}")

    order = [s for s in args.order.split(",") if s]
    names = discover(args.indir, order)
    if not names:
        sys.exit(f"No pll_1D_*.txt files found in {args.indir}")

    # ------------------------------------------------------------------
    # 1D confidence intervals from the profile scans
    # ------------------------------------------------------------------
    CI_LEVELS = {"68.3%": 1.00, "95%": 3.84}  # 1 dof thresholds
    print("\n1D confidence intervals (profile likelihood):")
    for name in names:
        x, y = load_1d(os.path.join(args.indir, f"pll_1D_{name}.txt"))
        best = x[np.nanargmin(y)]
        print(f"\n  {LABELS.get(name, name)}  (best fit = {best:.3g}):")
        for cl_name, thr in CI_LEVELS.items():
            cr = interval_1d(x, y, thr)
            if len(cr) >= 2:
                lo, hi = cr[0], cr[-1]
                print(f"    {cl_name}: [{lo:.4g}, {hi:.4g}]  "
                      f"(+{hi-best:.3g} / -{best-lo:.3g})")
            elif len(cr) == 1:
                print(f"    {cl_name}: one-sided bound at {cr[0]:.4g} "
                      f"(truncated by scan range)")
            else:
                print(f"    {cl_name}: not bracketed within scan range")

    n = len(names)
    print(f"Parameters ({n}): {', '.join(names)}")

    fig, axes = plt.subplots(n, n, figsize=(2.6 * n, 2.6 * n),
                             squeeze=False)

    for i in range(n):       # row
        for j in range(n):   # col
            ax = axes[i][j]

            # Upper triangle: blank
            if j > i:
                ax.axis("off")
                continue

            # Diagonal: 1D profile
            if i == j:
                name = names[i]
                x, y = load_1d(os.path.join(args.indir, f"pll_1D_{name}.txt"))
                ax.plot(x, y, color="black", lw=1.4)
                for lev in LEVELS_1D:
                    ax.axhline(lev, color="0.6", ls="--", lw=0.8)
                ax.set_ylim(0, max(6.0, np.nanpercentile(y, 95)))
                ax.set_xlim(x.min(), x.max())
                # mark the minimum
                xmin = x[np.argmin(y)]
                ax.axvline(xmin, color="tab:red", ls=":", lw=1.0)

            # Lower triangle: 2D contours
            else:
                xname = names[j]   # column -> x axis
                yname = names[i]   # row    -> y axis
                path, swapped = find_2d(args.indir, xname, yname)
                if path is None:
                    ax.text(0.5, 0.5, "no 2D", ha="center", va="center",
                            transform=ax.transAxes, color="0.6", fontsize=8)
                else:
                    X, Y, Z = load_2d(path)
                    if swapped:
                        # file stored (yname, xname): swap axes back
                        X, Y, Z = Y.T, X.T, Z.T
                    cs = ax.contour(X, Y, Z, levels=LEVELS_2D,
                                    colors=["tab:blue", "tab:cyan"],
                                    linewidths=1.2)
                    # best-fit point = grid minimum
                    iy, ix = np.unravel_index(np.nanargmin(Z), Z.shape)
                    ax.plot(X[iy, ix], Y[iy, ix], "r*", ms=8)

            # Axis labels only on edges
            if i == n - 1:
                ax.set_xlabel(LABELS.get(names[j], names[j]))
            else:
                ax.set_xticklabels([])
            if j == 0 and i != 0:
                ax.set_ylabel(LABELS.get(names[i], names[i]))
            elif j == 0 and i == 0:
                ax.set_ylabel(r"$2\Delta\ln L$")
            else:
                ax.set_yticklabels([])

            ax.tick_params(labelsize=8)

    # Legend (proxy artists)
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], color="tab:blue", lw=1.2, label=r"2D 68.3% ($2\Delta\ln L=2.30$)"),
        Line2D([0], [0], color="tab:cyan", lw=1.2, label=r"2D 95.4% ($2\Delta\ln L=6.18$)"),
        Line2D([0], [0], color="0.6", ls="--", lw=0.8, label=r"1D 68.3% / 95%"),
        Line2D([0], [0], color="tab:red", ls=":", lw=1.0, label="best fit"),
    ]
    fig.legend(handles=handles, loc="upper right", frameon=False, fontsize=10)

    fig.tight_layout()
    fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()