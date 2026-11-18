---
name: prism-statplot
description: "Prism-style group comparison for any continuous or ordinal readout - RNA expression, protein abundance (proteomics, western blot), qPCR, ELISA/concentration, cell proportions, scores. Chooses the statistical test from the design (2 groups vs 3+ groups, paired vs unpaired, parametric vs rank-based, all pairs vs against control), runs it with the matching multiple-comparison correction, and draws bar or dot + error bar + jittered replicates + stacked P-value brackets with a fixed group colour/marker style. Returns a methods sentence describing exactly what was done. The colour style is chosen in plain words from seven themes (bright, massier, npg, stallion, okabe, gastrulation, amp), e.g. \"use the massier style\", \"colour-blind safe\", \"换成 Nature 风格的配色\". Also gives one colour per cell type for single-cell data (UMAP, violin, bar): one hue per coarse lineage, shades for the fine cell types."
---

# prism-statplot: Prism-style statistical group comparison

One figure type for every "compare a readout across a few groups" question:
mean/median + error bar + every replicate as a point + significance brackets.
The test, the summary statistic, the y-axis floor and the colour roles are all
derived from **what the data is** and **how the experiment was designed**,
so the same rules apply to RNA, protein, qPCR or any other measurement.
The colours themselves come from a theme chosen in one line (section 6).

```python
# Inside Claude Science the functions are loaded automatically with the skill.
# In Claude Code or the Claude app: kernel.py is in this skill's folder; copy it
# to the working directory as prism_statplot.py, then
from prism_statplot import compare_groups, group_style, prism_group_plot, set_theme
# for 10 to 25 cell types also import celltype_palette (section 6)

set_theme("bright")       # colour theme for every plot; see section 6 for the choices
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
| `qpcr` | 2^-ΔΔCt fold change | **log2** (= -ΔΔCt), done automatically | parametric | bar | Relative expression (2^-ddCt) |
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

The design is passed to `compare_groups()` like this:

| Design | Arguments |
|---|---|
| every pair (default) | `comparisons="all"` |
| each group against a control | `comparisons="vs_control", control="Ctrl"` |
| paired / repeated | `paired=True` (values in the same subject order in every group) |
| clearly unequal variances | `equal_var=False` |
| force rank-based or parametric | `parametric=False` / `parametric=True` |

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

## 6. Colour themes and marker rules

A **theme** is a complete colour scheme: it colours every design below (two
groups, control + treatments, many groups, families) and the accent of
brackets, P values and titles. Choose it once, in one line:

```python
set_theme("massier")                      # every later plot uses it
style = group_style(order, theme="npg")   # or for one figure only
```

| Theme | Look | Source | Hues |
|---|---|---|---|
| `bright` (default) | pale pastel bars, saturated points, teal accent. The first two groups are the two-colour reference style: peach and salmon bars with coral points | this skill | 7 |
| `massier` | blue, soft yellow → orange, green, purple, brown, wine red | Massier et al. 2023, Nat Commun 14, 1438 (human white adipose tissue atlas) | 6 |
| `npg` | red, navy, teal, brown, cyan, salmon | ggsci "npg", modelled on figures in Nature Publishing Group journals | 6 |
| `stallion` | bold red, navy, green, purple, orange | ArchR "stallion" (Granja et al. 2021, Nat Genet), the cell-type colours of Kuppe et al. 2022, Nature 608, 766 | 7 |
| `okabe` | blue, vermilion, green, pink; colour-blind safe | Wu et al. 2021, Nat Genet 53, 1334 (breast cancer atlas), i.e. the Okabe-Ito set | 4 |
| `gastrulation` | soft, earthy: blue, red-orange, teal, lilac, rose, green | Pijuan-Sala et al. 2019, Nature 566, 490 (mouse gastrulation atlas) | 7 |
| `amp` | one light-to-dark ramp per lineage: blue, red, teal, magenta, brown | Zhang et al. 2019, Nat Immunol 20, 928 (rheumatoid arthritis synovium) | 5 |

All themes except `bright` use a near-black accent (`#222222`), as journal
figures do. `bright` keeps its teal accent (`#2D6869`). Error bars are black.

**Choosing the theme from what the user says.** Users pick a style in plain
words, in English or Chinese. Map their words to a theme and call `set_theme()`
once at the top of the script:

| The user says something like | Theme |
|---|---|
| "default", "pastel", "the original style", "peach and coral", "原来的风格", "两色风格", "淡彩" | `bright` |
| "massier", "the adipose atlas colours", "like the paper", "blue orange green purple", "文章那种配色", "脂肪图谱的配色" | `massier` |
| "Nature style", "NPG", "ggsci", "Nature 风格" | `npg` |
| "ArchR", "stallion", "bold", "vivid", "for a poster or slides", "鲜艳", "对比强" | `stallion` |
| "colour-blind safe", "Okabe-Ito", "色盲友好" | `okabe` |
| "soft", "earthy", "quiet", "柔和", "大地色", "素一点" | `gastrulation` |
| "one ramp per lineage", "AMP", "light to dark per cell class", "每个谱系一条色带" | `amp` |

- **No style mentioned:** use `bright`, and add one line to the reply that six
  other styles exist, with their names.
- **"Which styles are there?" / "有哪些配色":** list the seven themes with one
  line each. If the user already has data, offer to draw their figure in two or
  three of them side by side so they can choose.
- **"Switch to …" / "换成 …":** redraw the same figure with the new theme and
  change nothing else: same data, same test, same group order.
- **"Less saturated" / "饱和度低一点":** no theme has pale points, because points
  must stay visible on their bar. Use `gastrulation`, whose bars are the
  softest, and say that the points stay saturated on purpose.
- **A theme stays chosen.** Once the user picks a theme, use it for every later
  figure in the conversation until they pick another. Naming a style in any
  request, including one that only asks for a palette, counts as a new pick. For one figure in a
  different style, pass `theme=` to that call and leave the chosen theme alone.
- **A look that matches no theme** (for example "all grey" or a list of hex
  codes): use `binding={group: hex}` with the user's colours.
- State the theme used in the reply, by name.

**How a theme is built.** Each hue is three colours: light, dark and ink.

- A group's **bar** is a light shade between light and dark, its **points**
  are the dark colour. Bars have no outline.
- Where a source colour is too light to read as a point (cyan, orange), the
  point uses a deeper shade of the same hue and the source colour shows in the
  bars.
- When there are more groups than hues, the hues start a second round with a
  deeper bar and the ink colour for points, so no two groups look the same.
- The exact hex codes are in `THEMES` in `kernel.py`, with the source of each
  theme.

| Situation | `group_style(...)` | Result |
|---|---|---|
| groups with a control | `control="Ctrl"` | control takes the first hue, the others continue down the theme |
| ordered groups (dose, time) | `ordered=True` | multi-hue warm ramp yellow → orange → red → deep red (→ wine for 5+ groups); hue **and** lightness change together. The same ramp in every theme: a single-hue ramp of pastel bars differs by only ΔE ≈ 5 between neighbours |
| groups that already have a colour elsewhere (cell types) | `binding=palette` | the same hue as the UMAP / heatmap, bar = lightened tint |
| groups that form families (genotype × treatment, cell class × subtype, genotype × dose) | `families={group: family}` | one hue per family, members shaded light → dark in plotting order |
| none of the above | default | hues of the theme in order |

**Families: hue = family, shade = member.** For a two-level design the colour
shows the structure of the experiment.

- **Bars** step in lightness from the first to the last member (OKLab L 0.86 →
  0.68), so the last bar is clearly darker but never heavy, and every family
  steps by the same amount.
- **Points** run from the dark colour (first member) to the ink colour (last
  member). Bar and point darken together.
- **Limits.** As many families as the theme has hues (one more raises an
  error: merge groups or split the figure). At most 3 members per family; with
  4 the neighbouring shades differ by only ΔE ≈ 5.

```python
order = ["WT vehicle", "WT drug", "KO vehicle", "KO drug"]
style = group_style(order, families={"WT vehicle": "WT", "WT drug": "WT",
                                     "KO vehicle": "KO", "KO drug": "KO"})
```

**Measured contrast** (OKLab ΔE × 100; lowest value over all pairs of groups,
taking the better of bar and point; second number = simulated red-green colour
blindness):

| Theme | Point vs its bar | 4 groups | 2 × 2 families | 6 groups |
|---|---|---|---|---|
| `bright` | 18.6 | 10.1 / 2.8 | 19.6 / 16.5 | 10.1 / 2.8 |
| `massier` | 18.5 | 20.9 / 7.7 | 18.0 / 16.8 | 16.6 / 7.5 |
| `npg` | 19.0 | 15.5 / 4.8 | 18.0 / 13.8 | 8.3 / 4.8 |
| `stallion` | 23.4 | 19.5 / 4.9 | 18.4 / 12.1 | 11.0 / 4.9 |
| `okabe` | 18.3 | 16.4 / 7.6 | 18.8 / 17.0 | 16.4 / 7.6 |
| `gastrulation` | 19.7 | 16.7 / 6.1 | 18.3 / 17.3 | 11.4 / 3.1 |
| `amp` | 22.1 | 13.1 / 6.0 | 20.7 / 16.7 | 13.1 / 4.7 |

For readers with colour blindness use `okabe` or `massier`. The low value for
`bright` comes from its first two groups, which share the coral point colour by
design.

Every group also gets its own **marker shape** (o, s, ^, v, D, P …), so
groups stay distinguishable in grayscale, for colour-blind readers, and in the
two-group case of `bright` where both groups share the coral point colour.
Group order, theme and style are fixed once per figure series.

### Many cell types: `celltype_palette()`

`group_style()` is for the few groups of one plot. For the cell types of a
single-cell dataset, from a handful up to about 25, use `celltype_palette()`. It returns one colour per
cell type, `{cell type: hex}`: **one hue per lineage, shades within it**.

**Recommend two levels of annotation first.** The colours can only show a
structure that the annotation has. When a user has many cell types and a single
flat label, suggest annotating in two levels before choosing colours:

1. **Coarse**: cluster at low resolution and name five to seven major lineages
   from canonical markers (T/NK, myeloid, B, endothelial, stromal, epithelial).
2. **Fine**: re-cluster each lineage on its own and name the subtypes, so every
   fine cell type belongs to exactly one lineage. Keep both labels as two
   columns of `adata.obs`.

Then the hue shows the coarse level and the shade shows the fine level. Use the
coarse palette for overview figures and the fine palette for detailed ones:

```python
pairs = (adata.obs[["cell_type_coarse", "cell_type_fine"]]
         .drop_duplicates().sort_values(["cell_type_coarse", "cell_type_fine"]))
lineages = dict(zip(pairs["cell_type_fine"], pairs["cell_type_coarse"]))
fine_palette = celltype_palette(list(lineages), lineages)                   # shades
coarse_palette = celltype_palette(list(dict.fromkeys(lineages.values())))   # hues
```

Leaving out `lineages` colours the coarse level itself: each name takes the main
colour of its hue. List the coarse names in the order in which they first
appear among the fine cell types, so both palettes give a lineage the same hue.

The general form, with the mapping written by hand:

```python
lineages = {"CD4 naive T": "T/NK", "CD8 effector T": "T/NK", "NK": "T/NK",
            "Classical mono": "Myeloid", "Macrophage": "Myeloid",
            "Naive B": "B", "Plasma cell": "B", "Doublet": None}
palette = celltype_palette(list(lineages), lineages, theme="massier")

sc.pl.umap(adata, color="cell_type", palette=palette)       # scanpy
style = group_style(order, binding=palette)                 # the same colours in a bar plot
```

- Cell types are listed in legend order. Within a lineage the first one is the
  lightest and the last one the darkest.
- A lineage with one cell type takes the main colour of its hue. A cell type
  whose lineage is `None` is drawn light grey (`#D3D3D3`).
- **Shades differ in lightness and in hue.** Shades of one hue that differ only
  in lightness cannot be told apart as small points on a UMAP, and the two
  darkest ones both look black. So within a lineage:
  - lightness is evenly spaced between 0.86 and 0.42, at most 0.14 apart;
  - the hue shifts by up to 50° from the lightest to the darkest shade,
    centred on the lineage hue. Cool hues run from the cyan side to the violet
    side (blue: cyan → blue → indigo). Warm hues run from the yellow side to
    the red side. A ramp that already changes hue (yellow → orange → red in
    `massier`) is left as it is.
- **Limits.** As many lineages as the theme has hues (one more raises an
  error). Lowest ΔE between neighbouring shades by lineage size, over all
  themes: 2 or 3 cell types 14, 4 cell types 12 to 15, 5 cell types 9 to 12,
  6 cell types 7 to 9, 8 cell types 5 to 7. Keep a lineage to 5 cell types or
  fewer; split a larger one into two lineages (for example CD4 T and CD8 T/NK).
- For Seurat, copy the hex codes into a named vector for `DimPlot(cols = ...)`.

**Measured** on 20 cell types in 6 lineages (5/4/3/3/3/2), 190 pairs of
colours. A pair counts as hard to tell apart below ΔE 10 (below 8 under
simulated red-green colour blindness):

| Palette | Hard pairs | Lowest ΔE between neighbouring shades | Hard pairs, colour-blind |
|---|---|---|---|
| `bright` | 2 | 10.2 | 22 |
| `stallion` | 2 | 11.8 | 24 |
| `gastrulation` | 3 | 11.3 | 14 |
| `massier` | 5 | 11.1 | 23 |
| `npg` | 12 | 11.6 | 29 |
| `amp` (18 types, 5 lineages, 153 pairs) | 4 | 11.1 | 16 |
| `okabe` (15 types, 4 lineages, 105 pairs) | 0 | 11.2 | 6 |
| for comparison: ArchR stallion, 20 colours | 6 | – | 21 |
| for comparison: scanpy default, 20 colours | 12 | – | 32 |

Every theme except `npg` has fewer hard pairs than the two reference palettes.
No set of 20 colours keeps every pair apart, and none works for colour-blind
readers. Write the cell-type names on the clusters of a UMAP; do not rely on
the legend alone. ΔE understates the gain from the hue shift, so judge a
palette on a rendered UMAP as well.

## 7. Fixed geometry

- Bar width 0.6, **no outline** (no black frame, no coloured frame).
  Points size 38, filled, no rim. Error bar black, 1.2 pt, 4-pt caps,
  **drawn above the points** (`zorder` 4 vs 3).
- Points jittered within ±0.15; `seed` fixed per entity.
- Open box (top/right spines off). x labels `Group (n=x)`, regular weight,
  rotated 45° and right-anchored. y-axis title bold, in the accent colour of
  the theme.
- Brackets: thin line in the accent colour (1.0 pt), P value in regular weight,
  8.5 pt.
- The y-limits are computed from the data plus bracket headroom and never
  hard-coded.
- Save as PDF (vector) and PNG (300 dpi): `save_path` writes one file, so call
  `fig.savefig("name.png", dpi=300, bbox_inches="tight")` for the second format. Run a text-overlap check before
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
- More than about 8 groups: switch to a box/violin plot or heatmap. Their
  colours can still come from `celltype_palette()` (section 6).
- Differential expression across thousands of genes: use limma / DESeq2 and
  plot only selected genes here, with their model P values.
