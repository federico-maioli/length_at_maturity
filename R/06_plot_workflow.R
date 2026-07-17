library(here)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# 01 Define the workflow ----
# four phase-bands stacked vertically with extra spacing between them. Invisible left/right
# anchor points, aligned across bands, force every band (cluster) to the same width without
# changing node sizes. Edges are drawn as straight lines.

dot <- '
digraph workflow {
  graph [rankdir = TB, fontname = Helvetica, fontsize = 11,
         nodesep = 0.3, ranksep = 1.1, compound = true, newrank = true, splines = line]
  node  [shape = box, style = "rounded,filled", fontname = Helvetica,
         fontsize = 10, penwidth = 0, margin = "0.16,0.10"]
  edge  [color = grey45, arrowsize = 0.7]

  subgraph cluster_in {
    label = "Data source"; labeljust = l; fontname = "Helvetica-Bold";
    color = "#9DC3D6"; style = rounded;
    la_in [shape = point, width = 0.01, style = invis, group = "gleft"]
    ra_in [shape = point, width = 0.01, style = invis, group = "gright"]
    node [fillcolor = "#CDE7F0"]
    spawn  [label = "Spawning-season references\\n(GO-FISH, WKMAT, FishBase)"]
    survey [label = "ICES DATRAS survey data\\n(maturity records + haul locations)"]
  }

  subgraph cluster_proc {
    label = "Data standardization"; labeljust = l; fontname = "Helvetica-Bold";
    color = "#BDBDBD"; style = rounded;
    la_proc [shape = point, width = 0.01, style = invis, group = "gleft"]
    ra_proc [shape = point, width = 0.01, style = invis, group = "gright"]
    node [fillcolor = "#E6E6E6"]
    mclean [label = "data/final/maturity_clean.rds", shape = note, fillcolor = "#FDE3C4", penwidth = 1]
    season [label = "Filter to spawning /\\npre-spawning season;\\nlength-outlier checks"]
    harmon [label = "Harmonise maturity scales\\nto binary immature / mature;\\nassign ICES area & stock"]
  }

  subgraph cluster_mod {
    label = "Statistical analysis"; labeljust = l; fontname = "Helvetica-Bold";
    color = "#A9D3A9"; style = rounded;
    la_mod [shape = point, width = 0.01, style = invis, group = "gleft"]
    ra_mod [shape = point, width = 0.01, style = invis, group = "gright"]
    node [fillcolor = "#D5EFD5"]
    glm   [label = "Logistic GLM per stratum\\nL25 / L50 / L75\\nAUC, Tjur R2"]
    strat [label = "Stratify by scale & sex\\n(global / stock / area / period;\\nmin. 15 mature & immature)"]
  }

  subgraph cluster_out {
    label = "Technical validation"; labeljust = l; fontname = "Helvetica-Bold";
    color = "#D6A9A9"; style = rounded;
    la_out [shape = point, width = 0.01, style = invis, group = "gleft"]
    ra_out [shape = point, width = 0.01, style = invis, group = "gright"]
    valid [label = "FishBase plausibility check", fillcolor = "#F3D4D4"]
    final [label = "data/final/l50_estimates.csv", shape = note, fillcolor = "#FDE3C4", penwidth = 1]
    qc    [label = "Sanity filters\\n(positive, precise,\\nincreasing ogive, AUC >= 0.5)", fillcolor = "#E6E6E6"]
  }

  # each phase on one row, left and right anchors bracket the nodes
  {rank = same; la_in; survey; spawn; ra_in}
  {rank = same; la_proc; harmon; season; mclean; ra_proc}
  {rank = same; la_mod; strat; glm; ra_mod}
  {rank = same; la_out; qc; final; valid; ra_out}

  # strong pull of the first node to the left anchor packs the band to the left
  edge [constraint = false, style = invis, weight = 20]
  la_in   -> survey -> spawn
  la_proc -> harmon
  la_mod  -> strat
  la_out  -> qc

  # loose tie of the last node to the right anchor: the anchor floats out to the
  # widest band (equal width) while leaving the nodes packed at the left
  edge [constraint = false, style = invis, weight = 0]
  spawn  -> ra_in
  mclean -> ra_proc
  glm    -> ra_mod
  valid  -> ra_out

  # visible left-to-right flow within each band
  edge [constraint = false, style = solid]
  harmon -> season -> mclean
  strat  -> glm
  qc     -> final -> valid

  # both spines: align anchor columns across bands and stack the bands top-to-bottom
  edge [constraint = true, style = invis, weight = 100]
  la_in -> la_proc -> la_mod -> la_out
  ra_in -> ra_proc -> ra_mod -> ra_out

  # visible connectors, drawn as straight lines
  edge [constraint = false, style = solid, color = grey45]
  survey -> harmon
  spawn  -> season
  mclean -> strat
  glm    -> qc
}
'

# 02 Render and save ----

svg <- charToRaw(export_svg(grViz(dot)))

# vector pdf (best for publication) plus a high-resolution png (~300 dpi at full page width)
#rsvg_pdf(svg, here("outputs/supp/workflow.pdf"))
rsvg_png(svg, here("outputs/main/workflow.png"), width = 3200)
