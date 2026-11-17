---
name: prism-statplot
description: "Prism-style group comparison for any continuous or ordinal readout - RNA expression, protein abundance (proteomics, western blot), qPCR, ELISA/concentration, cell proportions, scores. Chooses the statistical test from the design (2 groups vs 3+ groups, paired vs unpaired, parametric vs rank-based, all pairs vs against control), runs it with the matching multiple-comparison correction, and draws bar or dot + error bar + jittered replicates + stacked P-value brackets with a fixed group colour/marker style. Returns a methods sentence describing exactly what was done."
---

# prism-statplot: Prism-style statistical group comparison

One figure type for every "compare a readout across a few groups" question:
mean/median + error bar + every replicate as a point + significance brackets.
The test, the summary statistic, the y-axis floor and the colours are all
derived from **what the data is** and **how the experiment was designed**,
so the same rules apply to RNA, protein, qPCR or any other measurement.

```python
# Inside Claude Science the functions are loaded automatically with the skill.
# Standalone: copy kernel.py next to your script as prism_statplot.py, then
from prism_statplot import compare_groups, group_style, prism_group_plot

res = compare_groups(values, order, data_type="protein", comparisons="all")
fig, ax = prism_group_plot(values, order,
                           style=group_style(order, control="Ctrl"),
                           brackets=res["brackets"], ylabel=res["ylabel"],
                           mode=res["mode"], summary=res["summary"],
                           title="ProteinX", save_path="ProteinX.pdf")
print(res["method"])      # paste into the figure legend / Methods
print(res["warnings"])    # read before trusting the P values
```

`compare_groups` is a convenience for the common designs below. When the data
were already modelled elsewhere (limma, DESeq2, a mixed model), pass those P
values as `brackets` directly — the model that accounts for the full design
always wins over a re-test on the plotted values.

## 1. What is one point? (decide this first)

**One point = one biological replicate** (one animal, one patient, one
independent culture). Technical replicates are averaged first. For
single-cell data, aggregate to one value per sample (pseudobulk mean or sum)
before testing — testing thousands of cells as if they were independent
inflates significance by orders of magnitude. `compare_groups` warns when a
group has n > 50 for this reason.

## 2. Data type → scale, test family, plot geometry

| `data_type` | Typical source | Test on | Default family | Geometry | y-label |
|---|---|---|---|---|---|
| `rna` | bulk RNA-seq TPM/CPM, scRNA pseudobulk | log2(x + 1), already transformed | parametric | **bar** (0 = not detected, a real zero) | Expression, log2(TPM + 1) |
| `protein` | LFQ / TMT / DIA proteomics | log2 intensity, already transformed | parametric | **dot** (0 on a log-intensity axis is meaningless; bars would exaggerate) | Protein abundance, log2 intensity |
| `qpcr` | 2^-ΔΔCt fold change | **log2** (= -ΔΔCt), done automatically | parametric | bar | Relative expression (2^-ΔΔCt) |
| `western` | densitometry / loading control | **log2**, automatic | parametric | bar | Relative band intensity |
| `concentration` | ELISA, cytokine, metabolite | **log2**, automatic | parametric | bar | Concentration (add units) |
| `proportion` | % cells, % area | as given | **rank-based** | bar | Cells (%) |
| `score` | histology / clinical ordinal | as given | **rank-based**, median shown | dot | Score |
| `generic` | anything else | as given | parametric | bar | — |

Why log-scale testing for ratios and concentrations: these readouts are
positive and right-skewed, and their biological effects are multiplicative
(2-fold up and 2-fold down should be symmetric). On log2 they are close to
normal, so a t-test or ANOVA is appropriate. The plot can still show the
linear values.

Rule for bars: draw a bar from zero **only when zero is a real zero**
(no expression, no protein, 0 %). Otherwise use `mode="dot"` with the axis
fitted to the data.

## 3. Design → test

| Groups | Design | Parametric (default for continuous) | Rank-based (proportions, scores, small skewed data) |
|---|---|---|---|
| 2 | unpaired | **Welch's t-test** (does not assume equal variance; default) | Mann-Whitney |
| 2 | paired (same subject in both) | paired t-test | Wilcoxon matched-pairs signed-rank |
| ≥3 | unpaired, all pairs | one-way ANOVA → **Tukey** | Kruskal-Wallis → **Dunn** (Holm-adjusted) |
| ≥3 | unpaired, each vs control | one-way ANOVA → **Dunnett** | Kruskal-Wallis → Dunn vs control |
| ≥3 | unpaired, clearly unequal variances (`equal_var=False`) | Welch's ANOVA → **Games-Howell** | — |
| ≥3 | paired / repeated | paired t-tests, Holm-adjusted | Friedman → Wilcoxon, Holm-adjusted |

Further rules:

- **Never run repeated t-tests across 3+ groups** without correction; the
  family-wise error grows with every comparison. The correction is part of the
  test choice above.
- **Choose all-pairs vs vs-control before looking at the data.** When the
  question is "does each treatment differ from control", Dunnett is the right
  test and has more power than Tukey.
- **Two-sided tests only**, alpha = 0.05.
- **Small n and rank tests.** With n = 3 vs 3, an exact Mann-Whitney test can
  never go below P = 0.10. With n < 5 per group, Dunn's test rests on a
  normal approximation and its P values are approximate. `compare_groups`
  warns in both cases. A parametric test on an appropriate scale (log2) is
  usually the better choice.
- **Two-factor designs** (genotype × treatment, cell type × time) are not
  one-way comparisons. Fit a two-way model (ANOVA, limma or a mixed model) and
  pass its P values in as `brackets`.

## 4. Summary statistic must match the test

| Test family | Centre | Error bar |
|---|---|---|
| parametric (t, ANOVA) | mean | **SEM** |
| rank-based | median | **interquartile range** |

`compare_groups` returns `summary` ("mean_sem" or "median_iqr"); pass it on so
the figure and the test never disagree. SD is preferable when the aim is to
show spread rather than precision; state which one you chose in the legend.

## 5. Which comparisons to annotate

- ≤ 3 comparisons: show all, including non-significant ones (exact P).
- More than 3: show only adjusted P < 0.05, and say so in the legend
  ("only significant comparisons are annotated"). `compare_groups` does this
  automatically and appends that sentence to `method`.
- P label: exact (`P = 0.0123`, `P < 0.0001`, `P > 0.9999`) by default.
  Use `label="stars"` (ns, \*, \*\*, \*\*\*, \*\*\*\*) only when the journal
  asks for it, and define the thresholds in the legend.
- Brackets are stacked automatically. The shortest spans go lowest and no two
  brackets share a level where they overlap in x. The vertical spacing is
  measured from the rendered label height, so labels never collide.

## 6. Colour and marker rules (bright style)

Pale, bright bars with saturated points; no outlines anywhere. Reference pair:

| Role | Bar | Points |
|---|---|---|
| group 1 / control | `#FAE2C5` peach | `#DA6F4E` coral, circle |
| group 2 | `#EFBAA0` salmon | `#DA6F4E` coral, square |
| groups 3-8 | `#CDE8F2` sky, `#D6EFD2` mint, `#E7DDF4` lavender, `#FFF0B8` lemon, `#F9D3E0` pink, `#D3F0EC` aqua | the saturated version of each hue |

Brackets, P values, the y-axis title and the plot title use one teal accent
(`#2D6869`). Error bars are black.

| Situation | `group_style(...)` | Result |
|---|---|---|
| groups with a control | `control="Ctrl"` | control takes the peach pair, the others continue down the palette |
| ordered groups (dose, time) | `ordered=True` | multi-hue warm ramp yellow → orange → red → deep red (→ wine for 5+ groups); hue **and** lightness change together. A single-hue ramp is not used: its pastel bars differ by only ΔE ≈ 5 and become indistinguishable |
| groups that already have a colour elsewhere (cell types) | `binding=palette` | the same hue as the UMAP / heatmap, bar = lightened tint |
| none of the above | default | palette in the order above |

Every group also gets its own **marker shape** (o, s, ^, v, D, P …), so
groups stay distinguishable in grayscale, for colour-blind readers, and in the
two-group case where both groups share the coral point colour.
Group order and style are fixed once per figure series.

## 7. Fixed geometry

- Bar width 0.6, **no outline** (no black frame, no coloured frame).
  Points size 38, filled, no rim. Error bar black, 1.2 pt, 4-pt caps,
  **drawn above the points** (`zorder` 4 vs 3).
- Points jittered within ±0.15; `seed` fixed per entity.
- Open box (top/right spines off). x labels `Group (n=x)`, regular weight,
  rotated 45° and right-anchored. y-axis title bold, teal.
- Brackets: thin teal line (1.0 pt), P value in regular weight, 8.5 pt.
- The y-limits are computed from the data plus bracket headroom and never
  hard-coded.
- Save as PDF (vector) and PNG (300 dpi). Run a text-overlap check before
  saving, and compare only tick labels that fall inside the visible axis range.

## 8. Legend / Methods sentence

`res["method"]` is a complete sentence. Example:

> Differences were assessed by one-way ANOVA followed by Dunnett's multiple
> comparisons test on log2-transformed values (two-sided, alpha = 0.05). Bars
> show mean +/- SEM; each point is one biological replicate (n = 4/4/4/4).

Add the software line once in Methods: statistics computed with SciPy
(`scipy.stats`, version as installed); plots with matplotlib.

## Not for

- Continuous x-axes (dose-response curves, time courses with many points):
  use a line plot or fit a curve.
- More than about 8 groups: switch to a box/violin plot or heatmap.
- Differential expression across thousands of genes: use limma / DESeq2 and
  plot only selected genes here, with their model P values.
