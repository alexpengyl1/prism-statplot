"""
prism-statplot - Prism-style group comparison: statistics chosen by design, summary matched to the
test, plot geometry matched to the data type. See SKILL.md for the decision rules.

compare_groups()   -> picks and runs the test, returns adjusted pairwise P values,
                      ready-made brackets, and a methods sentence
set_theme()        -> one line to switch the colour theme of every later plot
group_style()      -> per-group colours / marker shapes (control, ordered, bound, families)
celltype_palette() -> one colour per cell type for 10 to 25 cell types (hue = lineage)
prism_group_plot() -> draws the figure
"""

import math
import numpy as np

# What each data type is, how it is tested, and how it is drawn.
#   transform : applied before testing (None = test on the values as given)
#   parametric: default test family
#   mode      : "bar" only when zero on the axis is a real zero; otherwise "dot"
DATA_PRESETS = {'rna': {'transform': None,
         'parametric': True,
         'mode': 'bar',
         'ylabel': 'Expression, log2(TPM + 1)',
         'note': 'values already log2(TPM + 1) or log2(CPM + 1); 0 = not detected'},
 'protein': {'transform': None,
             'parametric': True,
             'mode': 'dot',
             'ylabel': 'Protein abundance, log2 intensity',
             'note': 'values already log2 intensity (LFQ / TMT / DIA); 0 is not a real zero, so no '
                     'bars'},
 'qpcr': {'transform': 'log2',
          'parametric': True,
          'mode': 'bar',
          'ylabel': 'Relative expression (2^-ddCt)',
          'note': 'values are 2^-ddCt fold changes; tested on log2 scale (= -ddCt)'},
 'concentration': {'transform': 'log2',
                   'parametric': True,
                   'mode': 'bar',
                   'ylabel': 'Concentration',
                   'note': 'positive linear readout (ELISA, cytokine, metabolite); tested on log2 '
                           'scale'},
 'western': {'transform': 'log2',
             'parametric': True,
             'mode': 'bar',
             'ylabel': 'Relative band intensity (normalised to loading control)',
             'note': 'densitometry ratio to loading control; tested on log2 scale'},
 'proportion': {'transform': None,
                'parametric': False,
                'mode': 'bar',
                'ylabel': 'Cells (%)',
                'note': 'bounded 0-100 %; rank-based tests by default'},
 'score': {'transform': None,
           'parametric': False,
           'mode': 'dot',
           'ylabel': 'Score',
           'note': 'ordinal score (histology, clinical); rank-based tests, median shown'},
 'generic': {'transform': None,
             'parametric': True,
             'mode': 'bar',
             'ylabel': '',
             'note': 'continuous readout tested as given'}}

# Bright pastel bar + saturated point. The first two pairs are the reference look
# (peach / salmon bars, coral points); further groups add equally bright hues.
# Same values as THEMES["bright"]["pairs"]; kept for scripts that import it.
DEFAULT_PALETTE = [
    ("#FAE2C5", "#DA6F4E"), ("#EFBAA0", "#DA6F4E"), ("#CDE8F2", "#3A93C3"),
    ("#D6EFD2", "#4FA85A"), ("#E7DDF4", "#8B62C2"), ("#FFF0B8", "#E0A800"),
    ("#F9D3E0", "#D9577F"), ("#D3F0EC", "#2BA69A"),
]
ORDERED_ANCHORS = ["#F5C84C", "#F2953A", "#E0603A", "#B8322E", "#7E1F2A"]
# Colour themes. One hue per group as (light, dark, ink):
#   bar   = a light shade between light and dark
#   point = dark; later members of a family step towards ink
# "pairs" lists ready-made (bar, point) pairs when a theme has them; otherwise each hue gives
# one mid-light bar and its dark point. "accent" colours brackets, P values and titles.
# Hues are ordered strongest contrast first. Where a source colour is too light to read as a
# point, the point uses a deeper shade of the same hue and the source colour shows in the bars.
THEMES = {
    # the original look: pale pastel bars, saturated points, teal accent;
    # the first two pairs are the two-colour reference style (peach / salmon bars, coral points)
    "bright": {
        "pairs": [("#FAE2C5", "#DA6F4E"), ("#EFBAA0", "#DA6F4E"), ("#CDE8F2", "#3A93C3"),
                  ("#D6EFD2", "#4FA85A"), ("#E7DDF4", "#8B62C2"), ("#FFF0B8", "#E0A800"),
                  ("#F9D3E0", "#D9577F"), ("#D3F0EC", "#2BA69A")],
        "ramps": [("#FAE2C5", "#DA6F4E", "#83432F"), ("#CDE8F2", "#3A93C3", "#235875"),
                  ("#D6EFD2", "#4FA85A", "#2F6536"), ("#E7DDF4", "#8B62C2", "#533B74"),
                  ("#FFF0B8", "#E0A800", "#866500"), ("#F9D3E0", "#D9577F", "#82344C"),
                  ("#D3F0EC", "#2BA69A", "#1A645C")],
        "accent": "#2D6869",
    },
    # Massier et al. 2023, Nat Commun 14, 1438 (human white adipose tissue atlas): blue, orange,
    # green and purple ramps from the paper's figure code, with the yellow end of the orange
    # ramp softened (the original #F3E55C is too saturated beside the other pale bars);
    # brown and wine red extend the set.
    "massier": {
        "pairs": None,
        "ramps": [("#C6DBEF", "#075A84", "#03304A"), ("#F7E08A", "#E8602D", "#A53A0C"),
                  ("#A6DBBB", "#359566", "#17513A"), ("#ECD9F1", "#967BCE", "#553A9A"),
                  ("#F0DBB9", "#815C0A", "#463204"), ("#FAD2D7", "#991B41", "#560D24")],
        "accent": "#222222",
    },
    # ggsci "npg": the palette modelled on figures in Nature Publishing Group journals
    # (nanxstats/ggsci, R/palettes.R). Six of its ten colours: red, navy, teal, brown, cyan, salmon.
    "npg": {
        "pairs": None,
        "ramps": [("#F8CDC6", "#E64B35", "#852C1F"), ("#D8DDE7", "#3C5488", "#23314F"),
                  ("#B8E4DD", "#00A087", "#005D4E"), ("#E5DFDA", "#7E6148", "#49382A"),
                  ("#CDECF3", "#3A8CA0", "#22515D"), ("#FCE3DB", "#AD6E5A", "#644034")],
        "accent": "#222222",
    },
    # ArchR "stallion" (Granja et al. 2021, Nat Genet 53, 403), the cell-type colours of
    # Kuppe et al. 2022, Nature 608, 766 (human myocardial infarction). Bold and saturated.
    # Seven of its twenty colours, in the palette's own order.
    "stallion": {
        "pairs": None,
        "ramps": [("#F3C0C2", "#D51F26", "#7C1216"), ("#D4D5E1", "#272E6A", "#171B3D"),
                  ("#C1DECA", "#208A42", "#135026"), ("#E7D4E9", "#89288F", "#4F1753"),
                  ("#FCDBC4", "#C06322", "#6F3914"), ("#CEE3E5", "#0C727C", "#074248"),
                  ("#EDD6E7", "#C06CAB", "#6F3F63")],
        "accent": "#222222",
    },
    # Wu et al. 2021, Nat Genet 53, 1334 (human breast cancer atlas): the eight major cell-type
    # colours, which are the Okabe-Ito colour-blind-safe set. The four that stay clearly apart
    # as point colours: blue, vermilion, green, pink.
    "okabe": {
        "pairs": None,
        "ramps": [("#CCE3F0", "#0072B2", "#004267"), ("#F3D2B8", "#D55E00", "#7C3700"),
                  ("#B8E4D8", "#009E73", "#005C43"), ("#F1D9E6", "#CC79A7", "#764661")],
        "accent": "#222222",
    },
    # Pijuan-Sala et al. 2019, Nature 566, 490 (mouse gastrulation atlas): seven of the 37
    # cell-type colours, picked for contrast. Soft and earthy, with a few strong accents.
    "gastrulation": {
        "pairs": None,
        "ramps": [("#CFDBEB", "#0F4A9C", "#092B5A"), ("#FBCDC1", "#EF4E22", "#8B2D14"),
                  ("#BDE2E0", "#139992", "#0B5955"), ("#DED7E8", "#8870AD", "#4F4164"),
                  ("#E4D4D3", "#9E6762", "#5C3C39"), ("#D4E7C9", "#65A83E", "#3B6124"),
                  ("#D7DCD3", "#354E23", "#1F2D14")],
        "accent": "#222222",
    },
    # Zhang et al. 2019, Nat Immunol 20, 928 (rheumatoid arthritis synovium): one light-to-dark
    # ramp per lineage (fibroblasts blue, B cells red, T cells teal, monocytes magenta), plus the
    # T-cell brown. The blue point colour is the midpoint of the paper's two darker blues.
    "amp": {
        "pairs": None,
        "ramps": [("#DEEBF7", "#3A6FA0", "#08306B"), ("#FCBBA1", "#CB181D", "#67000D"),
                  ("#C7EAE5", "#35978F", "#003C30"), ("#FDE0EF", "#AE017E", "#49006A"),
                  ("#E8DCCE", "#8C510A", "#512F06")],
        "accent": "#222222",
    },
}
THEME_STATE = {"name": "bright"}   # current theme; change it with set_theme()
ACCENT = "#2D6869"   # teal accent of the "bright" theme (kept for scripts that import it)
MARKER_SHAPES = ["o", "s", "^", "v", "D", "P", "X", "h"]


def format_pvalue(p, label="exact"):
    """'exact' -> Prism-style 'P = 0.0123'; 'stars' -> ns / * / ** / *** / ****."""
    if label == "stars":
        if p < 0.0001:
            return "****"
        if p < 0.001:
            return "***"
        if p < 0.01:
            return "**"
        if p < 0.05:
            return "*"
        return "ns"
    if p >= 0.9999:
        return "P > 0.9999"
    if p < 0.0001:
        return "P < 0.0001"
    return "P = {:.4f}".format(p)


def holm_adjust(pvals):
    """Holm step-down adjustment, order preserved."""
    p = np.asarray(pvals, dtype=float)
    order = np.argsort(p)
    m = len(p)
    adj = np.empty(m)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (m - rank) * p[idx])
        adj[idx] = min(running, 1.0)
    return adj


def assign_levels(pairs, x_positions):
    """Stack brackets so that no two brackets on the same level overlap in x.
    Shorter spans go lower, so nested comparisons read naturally."""
    spans = []
    for a, b in pairs:
        lo, hi = sorted((x_positions[a], x_positions[b]))
        spans.append((hi - lo, lo, hi, a, b))
    spans.sort()
    occupied = []
    levels = {}
    for _, lo, hi, a, b in spans:
        level = 0
        while any(lev == level and not (hi < l2 or lo > h2) for lev, l2, h2 in occupied):
            level += 1
        occupied.append((level, lo, hi))
        levels[(a, b)] = level
    return levels


def welch_anova(groups):
    """Welch's one-way ANOVA (unequal variances). Returns (F, df1, df2, P)."""
    from scipy import stats
    k = len(groups)
    n = np.array([len(g) for g in groups], dtype=float)
    m = np.array([np.mean(g) for g in groups])
    v = np.array([np.var(g, ddof=1) for g in groups])
    w = n / v
    grand = np.sum(w * m) / np.sum(w)
    a = np.sum(w * (m - grand) ** 2) / (k - 1)
    tmp = np.sum((1 - w / np.sum(w)) ** 2 / (n - 1))
    b = 1 + 2 * (k - 2) / (k ** 2 - 1) * tmp
    f = a / b
    df2 = (k ** 2 - 1) / (3 * tmp)
    return f, k - 1, df2, float(stats.f.sf(f, k - 1, df2))


def games_howell(groups, names):
    """Games-Howell pairwise test (unequal variances / unequal n)."""
    from scipy import stats
    k = len(groups)
    out = []
    for i in range(k):
        for j in range(i + 1, k):
            gi, gj = groups[i], groups[j]
            ni, nj = len(gi), len(gj)
            vi, vj = np.var(gi, ddof=1) / ni, np.var(gj, ddof=1) / nj
            q = abs(np.mean(gi) - np.mean(gj)) / math.sqrt((vi + vj) / 2.0)
            df = (vi + vj) ** 2 / (vi ** 2 / (ni - 1) + vj ** 2 / (nj - 1))
            out.append((names[i], names[j], float(stats.studentized_range.sf(q, k, df))))
    return out


def dunn_test(groups, names, control=None):
    """Dunn's rank-sum pairwise test after Kruskal-Wallis, tie-corrected, Holm-adjusted."""
    from scipy import stats
    pooled = np.concatenate(groups)
    ranks = stats.rankdata(pooled)
    N = len(pooled)
    _, counts = np.unique(pooled, return_counts=True)
    tie = np.sum(counts ** 3 - counts) / (12.0 * (N - 1))
    bounds = np.cumsum([0] + [len(g) for g in groups])
    mean_rank = [ranks[bounds[i]:bounds[i + 1]].mean() for i in range(len(groups))]
    pairs, raw = [], []
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            if control is not None and control not in (names[i], names[j]):
                continue
            ni, nj = len(groups[i]), len(groups[j])
            se = math.sqrt((N * (N + 1) / 12.0 - tie) * (1.0 / ni + 1.0 / nj))
            z = abs(mean_rank[i] - mean_rank[j]) / se
            pairs.append((names[i], names[j]))
            raw.append(2 * stats.norm.sf(z))
    adj = holm_adjust(raw)
    return [(a, b, float(p)) for (a, b), p in zip(pairs, adj)]


def min_rank_pvalue(n1, n2):
    """Smallest two-sided P an exact Mann-Whitney test can return for these n."""
    return 2.0 / math.comb(n1 + n2, n1)


def compare_groups(group_values, group_order, data_type="generic", paired=False,
                   parametric=None, comparisons="all", control=None, equal_var=None,
                   transform="preset", show=None, alpha=0.05):
    """Choose and run the right test for the design, return everything the plot needs.

    data_type   : key of DATA_PRESETS (sets transform, test family, plot mode)
    paired      : True when the same subject/sample appears in every group (same order)
    parametric  : override the preset test family
    comparisons : "all" (every pair) or "vs_control" (each group against `control`)
    equal_var   : None = Welch for 2 groups, Tukey/Dunnett for >=3; True / False to force
    show        : "all" or "significant"; default "all" when <= 3 pairs, else "significant"

    Returns dict(omnibus, pairs, brackets, summary, mode, ylabel, method, warnings)
    """
    from scipy import stats

    preset = DATA_PRESETS[data_type]
    if parametric is None:
        parametric = preset["parametric"]
    if transform == "preset":
        transform = preset["transform"]
    if comparisons == "vs_control" and control is None:
        control = group_order[0]

    warnings = []
    groups = []
    for g in group_order:
        vals = np.asarray(group_values[g], dtype=float)
        vals = vals[~np.isnan(vals)]
        if transform == "log2":
            if np.any(vals <= 0):
                raise ValueError("log2 transform needs positive values; group '%s' has <= 0" % g)
            vals = np.log2(vals)
        groups.append(vals)
    n = [len(g) for g in groups]
    k = len(groups)

    if min(n) < 3:
        warnings.append("a group has n < 3; P values are unreliable and SEM is not meaningful")
    if max(n) > 50:
        warnings.append("a group has n > 50; if points are cells rather than samples, aggregate "
                        "to one value per biological sample (pseudobulk) before testing")
    if paired and len(set(n)) > 1:
        raise ValueError("paired design needs equal n in every group (same subjects, same order)")
    if not parametric and not paired and k == 2:
        floor = min_rank_pvalue(n[0], n[1])
        if floor > alpha:
            warnings.append("Mann-Whitney cannot reach P < %.2f with n = %d/%d (smallest possible P = %.3f); "
                            "use parametric=True on a suitable scale, or add replicates" % (alpha, n[0], n[1], floor))
    if not parametric and not paired and k >= 3 and min(n) < 5:
        warnings.append("Dunn's test relies on a normal approximation that is unreliable below n = 5 per group; "
                        "treat its P values as approximate, or use parametric=True on a suitable scale")

    omnibus = None
    if k == 2:
        a, b = groups
        if paired and parametric:
            res, name = stats.ttest_rel(a, b), "paired t-test"
        elif paired:
            res, name = stats.wilcoxon(a, b), "Wilcoxon matched-pairs signed-rank test"
        elif parametric:
            ev = bool(equal_var) if equal_var is not None else False
            res = stats.ttest_ind(a, b, equal_var=ev)
            name = "unpaired t-test" if ev else "unpaired t-test with Welch's correction"
        else:
            res, name = stats.mannwhitneyu(a, b, alternative="two-sided"), "Mann-Whitney test"
        pairs = [(group_order[0], group_order[1], float(res.pvalue))]
        method = name
    elif paired:
        if parametric:
            omnibus = None
            raw, prs = [], []
            for i in range(k):
                for j in range(i + 1, k):
                    if control is not None and control not in (group_order[i], group_order[j]):
                        continue
                    prs.append((group_order[i], group_order[j]))
                    raw.append(stats.ttest_rel(groups[i], groups[j]).pvalue)
            method = "paired t-tests with Holm correction"
        else:
            fr = stats.friedmanchisquare(*groups)
            omnibus = dict(test="Friedman test", statistic=float(fr.statistic), p=float(fr.pvalue))
            raw, prs = [], []
            for i in range(k):
                for j in range(i + 1, k):
                    if control is not None and control not in (group_order[i], group_order[j]):
                        continue
                    prs.append((group_order[i], group_order[j]))
                    raw.append(stats.wilcoxon(groups[i], groups[j]).pvalue)
            method = "Friedman test followed by Wilcoxon signed-rank tests with Holm correction"
        pairs = [(a, b, float(p)) for (a, b), p in zip(prs, holm_adjust(raw))]
    elif parametric:
        use_welch = equal_var is False
        if use_welch:
            f, d1, d2, p = welch_anova(groups)
            omnibus = dict(test="Welch's ANOVA", statistic=float(f), p=p)
            pairs = games_howell(groups, group_order)
            if control is not None:
                pairs = [t for t in pairs if control in t[:2]]
            method = "Welch's ANOVA followed by Games-Howell multiple comparisons test"
        else:
            fa = stats.f_oneway(*groups)
            omnibus = dict(test="one-way ANOVA", statistic=float(fa.statistic), p=float(fa.pvalue))
            if comparisons == "vs_control":
                ci = group_order.index(control)
                others = [g for i, g in enumerate(groups) if i != ci]
                names = [g for i, g in enumerate(group_order) if i != ci]
                # Dunnett's P values come from a randomised integration; a fixed seed makes them
                # the same on every run (the argument is `rng` in SciPy >= 1.15, `random_state` before)
                import inspect
                seed_name = "rng" if "rng" in inspect.signature(stats.dunnett).parameters else "random_state"
                dn = stats.dunnett(*others, control=groups[ci], **{seed_name: 0})
                pairs = [(control, nm, float(p)) for nm, p in zip(names, dn.pvalue)]
                method = "one-way ANOVA followed by Dunnett's multiple comparisons test"
            else:
                tk = stats.tukey_hsd(*groups)
                pairs = [(group_order[i], group_order[j], float(tk.pvalue[i, j]))
                         for i in range(k) for j in range(i + 1, k)]
                method = "one-way ANOVA followed by Tukey's multiple comparisons test"
    else:
        kw = stats.kruskal(*groups)
        omnibus = dict(test="Kruskal-Wallis test", statistic=float(kw.statistic), p=float(kw.pvalue))
        pairs = dunn_test(groups, group_order, control=control if comparisons == "vs_control" else None)
        method = "Kruskal-Wallis test followed by Dunn's multiple comparisons test (Holm-adjusted)"

    if show is None:
        show = "all" if len(pairs) <= 3 else "significant"
    shown = [t for t in pairs if show == "all" or t[2] < alpha]
    x_positions = {g: i for i, g in enumerate(group_order)}
    levels = assign_levels([(a, b) for a, b, _ in shown], x_positions)
    brackets = [dict(group_a=a, group_b=b, pvalue=p, level=levels[(a, b)]) for a, b, p in shown]

    scale = " on log2-transformed values" if transform == "log2" else ""
    summary = "mean_sem" if parametric else "median_iqr"
    shown_txt = "mean +/- SEM" if parametric else "median with interquartile range"
    method_sentence = ("Differences were assessed by {}{} (two-sided, alpha = {}). Bars show {}; "
                       "each point is one biological replicate (n = {}).").format(
        method, scale, alpha, shown_txt, "/".join(str(x) for x in n))
    if show == "significant" and len(shown) < len(pairs):
        method_sentence += " Only comparisons with adjusted P < {} are annotated.".format(alpha)

    return dict(omnibus=omnibus, pairs=pairs, brackets=brackets, summary=summary,
                mode=preset["mode"], ylabel=preset["ylabel"], method=method_sentence,
                n=dict(zip(group_order, n)), warnings=warnings)


def lighten(hex_color, amount=0.75):
    """Mix a colour with white; amount 0 = unchanged, 1 = white."""
    h = hex_color.lstrip("#")
    rgb = [int(h[i:i + 2], 16) for i in (0, 2, 4)]
    mixed = [int(round(c + (255 - c) * amount)) for c in rgb]
    return "#{:02X}{:02X}{:02X}".format(*mixed)


def mix(hex_a, hex_b, t):
    """Blend two colours in RGB; t = 0 gives hex_a, t = 1 gives hex_b."""
    a = [int(hex_a[k:k + 2], 16) for k in (1, 3, 5)]
    b = [int(hex_b[k:k + 2], 16) for k in (1, 3, 5)]
    return "#{:02X}{:02X}{:02X}".format(*[int(round(x + (y - x) * t)) for x, y in zip(a, b)])


def oklab(hex_color):
    """OKLab coordinates (L, a, b) of a colour; L runs from 0 = black to 1 = white."""
    srgb = [int(hex_color[k:k + 2], 16) / 255 for k in (1, 3, 5)]
    r, g, b = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in srgb]
    l = (0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b) ** (1 / 3)
    m = (0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b) ** (1 / 3)
    s = (0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b) ** (1 / 3)
    return (0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
            1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
            0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s)


def lightness(hex_color):
    """OKLab lightness of a colour: 0 = black, 1 = white."""
    return oklab(hex_color)[0]


def oklch(hex_color):
    """(lightness, chroma, hue angle in degrees) of a colour."""
    L, a, b = oklab(hex_color)
    return L, math.hypot(a, b), math.degrees(math.atan2(b, a)) % 360


def from_oklch(L, chroma, hue):
    """Hex colour with this lightness, chroma and hue; chroma is reduced until the colour fits in sRGB."""
    def linear_rgb(c):
        a, b = c * math.cos(math.radians(hue)), c * math.sin(math.radians(hue))
        l = (L + 0.3963377774 * a + 0.2158037573 * b) ** 3
        m = (L - 0.1055613458 * a - 0.0638541728 * b) ** 3
        s = (L - 0.0894841775 * a - 1.2914855480 * b) ** 3
        return (4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
                -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
                -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s)

    def fits(c):
        return all(-1e-4 <= x <= 1 + 1e-4 for x in linear_rgb(c))

    if not fits(chroma):
        lo, hi = 0.0, chroma
        for _ in range(20):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if fits(mid) else (lo, mid)
        chroma = lo
    rgb = [min(max(x, 0.0), 1.0) for x in linear_rgb(chroma)]
    rgb = [12.92 * x if x <= 0.0031308 else 1.055 * x ** (1 / 2.4) - 0.055 for x in rgb]
    return "#{:02X}{:02X}{:02X}".format(*[int(round(255 * x)) for x in rgb])


def shade_at(light, dark, target):
    """Colour between `light` and `dark` whose lightness is `target` (clamped to the two ends)."""
    lo, hi = 0.0, 1.0
    for _ in range(20):
        mid = (lo + hi) / 2
        if lightness(mix(light, dark, mid)) > target:
            lo = mid
        else:
            hi = mid
    return mix(light, dark, (lo + hi) / 2)


def set_theme(name):
    """Choose the colour theme used by every later group_style() / prism_group_plot() call."""
    if name not in THEMES:
        raise KeyError("unknown theme '{}'; choose one of: {}".format(name, ", ".join(THEMES)))
    THEME_STATE["name"] = name


def group_style(group_order, control=None, ordered=False, binding=None, families=None, theme=None):
    """Per-group dict(bar, marker, marker_shape, edge).

    theme   : name of a colour theme in THEMES for this call; default = the one set by set_theme()

    control : this group takes the first hue of the theme; the others continue down the theme
    ordered : groups are a dose / time series -> one hue ramp, light to dark
    binding : {group: hex} already fixed elsewhere (e.g. cell-type colours) -> reused
    families: {group: family} for two-level designs (genotype x treatment, cell class x subtype)
              -> one hue per family, members shaded light to dark in plotting order.
              Up to 3 members per family keep neighbouring shades clearly apart.
    """
    colours = THEMES[theme or THEME_STATE["name"]]
    ramps = colours["ramps"]
    # ready-made (bar, point) pairs, or per hue one mid-light bar with the dark point;
    # when the groups outnumber the hues, a second round uses a deeper bar with the ink point
    pairs = colours["pairs"] or ([(mix(light, dark, 0.25), dark) for light, dark, ink in ramps] +
                                 [(shade_at(light, dark, 0.70), ink) for light, dark, ink in ramps])
    style = {}
    if families is not None:
        names = list(dict.fromkeys(families[g] for g in group_order))   # hue order = first appearance
        if len(names) > len(ramps):
            raise ValueError("families: this theme has {} hues; merge some families or split the figure".format(
                len(ramps)))
        for f, name in enumerate(names):
            light, dark, ink = ramps[f]
            members = [g for g in group_order if families[g] == name]
            m = len(members)
            for j, g in enumerate(members):
                if m == 1:      # a single member: one mid-light bar and the dark point
                    bar, point = mix(light, dark, 0.25), dark
                else:           # bars step from lightness 0.86 down to 0.68, points from dark to ink
                    bar = shade_at(light, dark, 0.86 - 0.18 * j / (m - 1))
                    point = mix(dark, ink, j / (m - 1))
                style[g] = dict(bar=bar, marker=point,
                                marker_shape=MARKER_SHAPES[group_order.index(g) % len(MARKER_SHAPES)],
                                edge=None)
        return style
    if binding is not None:
        for i, g in enumerate(group_order):
            style[g] = dict(bar=lighten(binding[g], 0.72), marker=binding[g],
                            marker_shape=MARKER_SHAPES[i % len(MARKER_SHAPES)], edge=None)
        return style
    if ordered:
        # Multi-hue sequential ramp (yellow -> orange -> red -> deep red -> wine): hue and
        # lightness change together, so neighbouring doses stay distinguishable.
        # A single-hue ramp collapses to near-identical pastel bars (adjacent dE ~ 5).
        anchors = ORDERED_ANCHORS[:4] if len(group_order) <= 4 else ORDERED_ANCHORS
        rgb = [np.array([int(a[k:k + 2], 16) for k in (1, 3, 5)], dtype=float) for a in anchors]
        pos = np.linspace(0, len(rgb) - 1, len(group_order))
        for i, g in enumerate(group_order):
            seg = min(int(pos[i]), len(rgb) - 2)
            frac = pos[i] - seg
            c = rgb[seg] * (1 - frac) + rgb[seg + 1] * frac
            mk = "#{:02X}{:02X}{:02X}".format(*[int(round(x)) for x in c])
            style[g] = dict(bar=lighten(mk, 0.50), marker=mk,
                            marker_shape=MARKER_SHAPES[i % len(MARKER_SHAPES)], edge=None)
        return style
    j = 0
    for i, g in enumerate(group_order):
        if g == control:
            bar, mk = pairs[0]
            style[g] = dict(bar=bar, marker=mk, marker_shape="o", edge=None)
            continue
        start = 1 if control in group_order else 0
        bar, mk = pairs[(j + start) % len(pairs)]
        style[g] = dict(bar=bar, marker=mk, marker_shape=MARKER_SHAPES[(j + start) % len(MARKER_SHAPES)],
                        edge=None)
        j += 1
    return style


def celltype_palette(celltypes, lineages=None, theme=None, other="#D3D3D3"):
    """{cell type: hex} for many cell types (about 10 to 25): one hue per lineage, shades within it.

    celltypes: cell types in legend order; within a lineage the first is the lightest
    lineages : {cell type: lineage}, e.g. {"CD4 T": "T/NK", "CD8 T": "T/NK", "B": "B"};
               a cell type whose lineage is None takes the grey `other` colour.
               Leave it out to colour the coarse level itself: each name then takes the main
               colour of its own hue, e.g. celltype_palette(["T/NK", "Myeloid", "B"]).
               Hues are given in order of first appearance, so list the coarse names in the
               order in which they first appear among the fine cell types.
    theme    : name of a colour theme in THEMES for this call; default = the one set by set_theme()

    The dict goes straight into scanpy (palette=), seaborn (palette=) and
    group_style(order, binding=...), so UMAP, violin and bar plots share one colour per cell type.

    Within a lineage the shades differ in lightness AND in hue, because shades of one hue that
    differ only in lightness cannot be told apart as small points on a UMAP:
      lightness: evenly spaced between 0.86 and 0.42, at most 0.14 apart (no near-black shades)
      hue      : shifts by up to 50 degrees from the lightest to the darkest shade, centred on the
                 lineage hue. Cool hues run from the cyan side to the violet side (blue: cyan ->
                 blue -> indigo); warm hues from the yellow side to the red side.
                 A ramp that already changes hue (yellow -> orange -> red) is left as it is.
    Up to 5 cell types per lineage stay clearly apart.
    """
    if lineages is None:      # coarse level: every name is its own lineage; None (unassigned) is skipped
        celltypes = [c for c in celltypes if c is not None]
        lineages = {c: c for c in celltypes}
    ramps = THEMES[theme or THEME_STATE["name"]]["ramps"]
    names = list(dict.fromkeys(lineages[c] for c in celltypes if lineages[c] is not None))
    if len(names) > len(ramps):
        raise ValueError("celltype_palette: this theme has {} hues but there are {} lineages; "
                         "merge some lineages or choose a theme with more hues".format(len(ramps), len(names)))
    palette = {c: other for c in celltypes}
    for f, name in enumerate(names):
        light, dark, ink = ramps[f]
        members = [c for c in celltypes if lineages[c] == name]
        m = len(members)
        if m == 1:      # a single cell type takes the main colour of its hue
            palette[members[0]] = dark
            continue
        top, bottom = min(0.86, lightness(light)), max(0.42, lightness(ink))
        step = min((top - bottom) / (m - 1), 0.14)
        first = (top + bottom) / 2 + step * (m - 1) / 2      # shades are centred in the range
        dark_chroma, dark_hue = oklch(dark)[1:]
        own_change = (oklch(ink)[2] - oklch(light)[2] + 180) % 360 - 180    # hue change the ramp already has
        sweep = 0 if abs(own_change) > 25 else 50 * min(1, (m - 1) / 4)     # smaller for 2 to 4 cell types
        direction = 1 if 110 <= dark_hue <= 330 else -1                     # cool hues / warm hues
        for j, c in enumerate(members):
            target = first - step * j
            shade = shade_at(light, dark, target) if target >= lightness(dark) else shade_at(dark, ink, target)
            if sweep:
                position = j / (m - 1) - 0.5                 # -0.5 lightest ... +0.5 darkest
                shade = from_oklch(target, max(oklch(shade)[1], dark_chroma),
                                   (dark_hue + direction * position * sweep) % 360)
            palette[c] = shade
    return palette


def prism_group_plot(group_values, group_order, style=None, brackets=None, ylabel="",
                     title="", save_path=None, seed=0, figsize=None, title_color=None,
                     mode="bar", summary="mean_sem", label="exact", show_n=True, theme=None):
    """Bar (or dot) + error bar + jittered replicates + stacked significance brackets.

    theme   : colour theme for the accent (brackets, P values, titles) and the default style

    mode    : "bar" draws a filled bar from 0 (only when 0 is a real zero);
              "dot" draws a short mean/median line, y-axis fitted to the data
    summary : "mean_sem" or "median_iqr" - must match the test (use compare_groups output)
    label   : "exact" (P = 0.0123) or "stars"
    show_n  : append (n = x) under each group label
    """
    import matplotlib.pyplot as plt

    if title_color is None:
        title_color = THEMES[theme or THEME_STATE["name"]]["accent"]
    if style is None:
        style = group_style(group_order, theme=theme)
    brackets = brackets or []
    missing = [g for g in group_order if g not in group_values]
    if missing:
        raise KeyError("group_values missing entries for: {}".format(missing))
    if figsize is None:
        figsize = (1.4 + 1.0 * len(group_order), 4.6)

    rng = np.random.default_rng(seed)
    xs = {g: i for i, g in enumerate(group_order)}
    fig, ax = plt.subplots(figsize=figsize)

    values = {}
    for g in group_order:
        v = np.asarray(group_values[g], dtype=float)
        v = v[~np.isnan(v)]
        values[g] = v
        sty, x0 = style[g], xs[g]
        if summary == "median_iqr":
            centre = float(np.median(v))
            q1, q3 = np.percentile(v, [25, 75])
            err = [[centre - q1], [q3 - centre]]
        else:
            centre = float(np.mean(v))
            err = [[np.std(v, ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0]] * 2
        if mode == "bar":
            edge = sty.get("edge")
            ax.bar(x0, centre, width=0.6, color=sty["bar"], edgecolor=edge if edge else "none",
                   linewidth=1.3 if edge else 0, zorder=1)
        else:
            ax.plot([x0 - 0.25, x0 + 0.25], [centre, centre], color="black", linewidth=1.8, zorder=4)
        ax.errorbar(x0, centre, yerr=err, fmt="none", ecolor="black", elinewidth=1.2, capsize=4, zorder=4)
        jitter = rng.uniform(-0.15, 0.15, size=len(v))
        ax.scatter(x0 + jitter, v, facecolor=sty["marker"], edgecolor=sty["marker"],
                   marker=sty["marker_shape"], s=38, linewidths=0, zorder=3)

    labels = ["{} (n={})".format(g, len(values[g])) if show_n else g for g in group_order]
    ax.set_xticks(list(xs.values()))
    ax.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor", fontsize=9)
    ax.set_ylabel(ylabel, fontsize=10.5, fontweight="bold", color=title_color)
    if title:
        ax.set_title(title, fontsize=13, color=title_color, fontweight="bold")

    allv = np.concatenate([values[g] for g in group_order])
    lo, hi = float(allv.min()), float(allv.max())
    yr = hi - lo if hi > lo else (abs(hi) or 1.0)
    floor = min(0.0, lo - yr * 0.08) if mode == "bar" else lo - yr * 0.12
    tick = yr * 0.03
    n_levels = max((b["level"] for b in brackets), default=-1) + 1

    def draw(step, base):
        arts = []
        for b in brackets:
            x1, x2 = xs[b["group_a"]], xs[b["group_b"]]
            y = base + b["level"] * step
            line, = ax.plot([x1, x1, x2, x2], [y, y + tick, y + tick, y], color=title_color, linewidth=1.0)
            txt = ax.text((x1 + x2) / 2, y + tick * 1.6, format_pvalue(b["pvalue"], label),
                          ha="center", va="bottom", fontsize=8.5, color=title_color)
            arts.extend([line, txt])
        return arts

    if brackets:
        step, base = yr * 0.16, hi + yr * 0.10
        ax.set_ylim(floor, hi + yr * (0.18 + 0.22 * n_levels))
        arts = draw(step, base)
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        inv = ax.transData.inverted()
        th = 0.0
        for a in arts:
            if hasattr(a, "get_text"):
                bb = a.get_window_extent(r)
                th = max(th, abs(inv.transform((0, bb.y1))[1] - inv.transform((0, bb.y0))[1]))
        need = th * 1.8 + tick * 1.6
        if need > step:
            for a in arts:
                a.remove()
            step = need
            ax.set_ylim(floor, hi + step * n_levels + yr * 0.18)
            draw(step, base)
    else:
        ax.set_ylim(floor, hi + yr * 0.18)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", labelsize=10)
    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig, ax


def default_group_style(group_order):
    """Backward-compatible name from v1; same as group_style(group_order)."""
    return group_style(group_order)
