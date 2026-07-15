"""
Shared utilities for curve fitting and plot styling.
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d


# ─── Hill equation ───────────────────────────────────────────────────────────

def hill4(x, baseline, df_max, ec50, n):
    """4-parameter Hill / logistic equation."""
    return baseline + df_max / (1.0 + (ec50 / x) ** n)


NOISE_THRESHOLD = 0.12   # dF/F0 range below which data is treated as noise
N_BOUNDARY = 9.5         # Hill n at or above this is considered a boundary hit


def fit_hill(conc, signal, sem=None, n_starts=6):
    """
    Fit 4-parameter Hill equation to dose-response data.

    Returns dict with keys: baseline, df_max, ec50, n, r2, success, message,
                             is_noise (bool), is_step_fn (bool)
    """
    x = np.asarray(conc, dtype=float)
    y = np.asarray(signal, dtype=float)

    signal_range = y.max() - y.min()
    # Flag datasets where the signal range is too small to be meaningful
    is_noise = signal_range < NOISE_THRESHOLD

    # Flag step-function data: ≥6 of first 7 points within noise range
    # but the last point is notably higher (instrument-limited/single-point response)
    is_step_fn = False
    if len(y) >= 4:
        # Check if ≥80% of all-but-last points have range < 0.15
        body_range = y[:-1].max() - y[:-1].min()
        last_jump  = y[-1] - y[:-1].mean()
        if body_range < NOISE_THRESHOLD and last_jump > 0.2:
            is_step_fn = True

    # Sensible initial guesses
    ymin, ymax = y.min(), y.max()
    xmid = np.median(x)

    starts = []
    for ec50_0 in [xmid, x[len(x)//2], x.max()/2, x.max()*2]:
        for n_0 in [1.0, 2.0]:
            starts.append([ymin, ymax - ymin, ec50_0, n_0])

    bounds = (
        [-np.inf, 0,      1e-4,  0.1],
        [ np.inf, np.inf, 1e6,  10.0],
    )

    best = None
    best_r2 = -np.inf

    for p0 in starts:
        try:
            popt, _ = curve_fit(
                hill4, x, y, p0=p0, bounds=bounds,
                maxfev=10000, method='trf'
            )
            y_pred = hill4(x, *popt)
            ss_res = np.sum((y - y_pred) ** 2)
            ss_tot = np.sum((y - y.mean()) ** 2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
            if r2 > best_r2:
                best_r2 = r2
                best = popt
        except Exception:
            pass

    if best is None:
        return dict(baseline=np.nan, df_max=np.nan, ec50=np.nan, n=np.nan,
                    r2=np.nan, success=False, message='fit failed')

    baseline, df_max, ec50, n = best

    # Determine quality flags and message
    flags = []
    if is_noise:
        flags.append('noise: signal range < 0.15 dF/F0')
    if is_step_fn:
        flags.append('step-function: EC50 unreliable')
    if n >= N_BOUNDARY:
        flags.append(f'Hill n at boundary ({n:.1f})')
    msg = '; '.join(flags) if flags else 'ok'

    return dict(baseline=baseline, df_max=df_max, ec50=ec50, n=n,
                r2=best_r2, success=True, message=msg,
                is_noise=is_noise, is_step_fn=is_step_fn)


def compute_lod(conc, fit_params, sem_values):
    """
    LOD = concentration where fitted signal = baseline + 3 * mean(|SEM|).
    Interpolated from a dense grid of the fit curve.
    Returns LOD in µM, or np.nan if not found.
    """
    baseline = fit_params['baseline']
    df_max   = fit_params['df_max']
    ec50     = fit_params['ec50']
    n        = fit_params['n']

    sem_arr = np.abs(np.asarray(sem_values, dtype=float))
    sem_mean = sem_arr[sem_arr > 0].mean() if np.any(sem_arr > 0) else 0.0

    threshold = baseline + 3 * sem_mean

    # Dense grid from 1e-4 to max(conc)*10
    x_dense = np.logspace(np.log10(1e-4), np.log10(max(conc) * 10), 5000)
    y_dense = hill4(x_dense, baseline, df_max, ec50, n)

    # Find first crossing
    above = y_dense >= threshold
    if not np.any(above):
        return np.nan
    idx = np.argmax(above)
    if idx == 0:
        return x_dense[0]
    # Linear interpolation between points straddling threshold
    x0, x1 = x_dense[idx - 1], x_dense[idx]
    y0, y1 = y_dense[idx - 1], y_dense[idx]
    lod = x0 + (threshold - y0) * (x1 - x0) / (y1 - y0)
    return lod


# ─── Plot styling ─────────────────────────────────────────────────────────────

CATEGORY_COLORS = {
    'Environmental':    '#2ca02c',
    'Hormones/Steroids':'#9467bd',
    'Pharmaceuticals':  '#1f77b4',
    'Stimulants':       '#ff7f0e',
    'Metabolites':      '#8c564b',
    'Neurotransmitters':'#e377c2',
    'Amino Acids':      '#bcbd22',
    'Other':            '#7f7f7f',
}

def publication_style():
    """Apply publication-quality matplotlib style (Nature/ACS aesthetic)."""
    matplotlib.rcParams.update({
        'font.family':        'sans-serif',
        'font.sans-serif':    ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size':          8,
        'axes.labelsize':     9,
        'axes.titlesize':     9,
        'xtick.labelsize':    8,
        'ytick.labelsize':    8,
        'legend.fontsize':    7,
        'axes.linewidth':     0.8,
        'axes.spines.top':    False,
        'axes.spines.right':  False,
        'xtick.major.width':  0.8,
        'ytick.major.width':  0.8,
        'xtick.minor.width':  0.5,
        'ytick.minor.width':  0.5,
        'lines.linewidth':    1.2,
        'lines.markersize':   4,
        'figure.dpi':         150,
        'savefig.dpi':        300,
        'savefig.bbox':       'tight',
        'pdf.fonttype':       42,   # embed fonts
        'ps.fonttype':        42,
        'axes.unicode_minus': True,
    })
    # Use mathtext for subscripts so they render in any font
    matplotlib.rcParams['text.usetex'] = False


def save_figure(fig, path_stem):
    """Save figure as both PNG and PDF."""
    fig.savefig(path_stem + '.png', dpi=300, bbox_inches='tight')
    fig.savefig(path_stem + '.pdf', dpi=300, bbox_inches='tight')
    print(f"  Saved {path_stem}.png / .pdf")


def get_sensor_colors(n):
    """Return n color-blind-friendly colors for sensor curves."""
    cmap = matplotlib.colormaps.get_cmap('tab10').resampled(max(n, 10))
    return [cmap(i) for i in range(n)]


# ─── Chemical category map ────────────────────────────────────────────────────

CATEGORY_MAP = {
    # Environmental contaminants
    'DEHP':           'Environmental',
    'MEHP':           'Environmental',
    'Phthalic acid':  'Environmental',
    'Fipronil':       'Environmental',
    'Atrazine':       'Environmental',
    'Glyphosate':     'Environmental',
    'DEET':           'Environmental',
    # Hormones / steroids
    'Estradiol':         'Hormones/Steroids',
    'Progesterone':      'Hormones/Steroids',
    'DHEA':              'Hormones/Steroids',
    'Corticosterone':    'Hormones/Steroids',
    'L-Thyroxine':       'Hormones/Steroids',
    'T3':                'Hormones/Steroids',
    # Pharmaceuticals
    'Ciprofloxacin':     'Pharmaceuticals',
    'Atenolol':          'Pharmaceuticals',
    'Cyclosporin A':     'Pharmaceuticals',
    'Everolimus':        'Pharmaceuticals',
    'Carbamazepine':     'Pharmaceuticals',
    'Phenytoin':         'Pharmaceuticals',
    'Lamotrigine':       'Pharmaceuticals',
    'Valproic acid':     'Pharmaceuticals',
    'Digoxin':           'Pharmaceuticals',
    'Colchicine':        'Pharmaceuticals',
    'Glyburide':         'Pharmaceuticals',
    'Rosuvastatin':      'Pharmaceuticals',
    'Atorvastatin':      'Pharmaceuticals',
    # Stimulants / nicotinic
    'Nicotine':          'Stimulants',
    'Caffeine':          'Stimulants',
    'Theophylline':      'Stimulants',
    'Theobromine':       'Stimulants',
    'Ritalinic acid':    'Stimulants',
    # Metabolites / vitamins
    'Thiamine':          'Metabolites',
    'Ergothioneine':     'Metabolites',
    'Carnitine':         'Metabolites',
    'Beta-Hydroxybutyrate': 'Metabolites',
    'Glutathione':       'Metabolites',
    'Folic acid':        'Metabolites',
    'Nicotinamide':      'Metabolites',
    'Taurine':           'Metabolites',
    'Creatinine':        'Metabolites',
    'Uric acid':         'Metabolites',
    'Citrulline':        'Metabolites',
    'Melatonin':         'Metabolites',
    'Aspartame':         'Metabolites',
    'D-Neopterin':       'Metabolites',
    'Bilirubin':         'Metabolites',
    # Neurotransmitters / amines
    'Dopamine':          'Neurotransmitters',
    'L-Dopa':            'Neurotransmitters',
    'Carbidopa':         'Neurotransmitters',
    'Norepinephrine':    'Neurotransmitters',
    'Epinephrine':       'Neurotransmitters',
    'Tyramine':          'Neurotransmitters',
    'Tryptamine':        'Neurotransmitters',
    'Histamine':         'Neurotransmitters',
    'Betahistine':       'Neurotransmitters',
    # Amino acids
    'Carnosine':         'Amino Acids',
    'Theanine':          'Amino Acids',
    '3-Methyl-L-histidine': 'Amino Acids',
    # Other
    'Glucosamine':       'Other',
    'Sodium lactate':    'Other',
    'MST':               'Other',
    'Milrinone':         'Other',
    'Acetaminophen':     'Other',
}
