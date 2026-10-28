# Angstel/Vecht bog-oak chronology: code and data

Arts, M. and Bazelmans, J., *Decadal to multidecadal variability in a nearly
1.5-millennia-long bog-oak chronology from the Angstel/Vecht delta, the Netherlands
(1283 BCE–156 CE)*. Submitted to *The Holocene*.

## How to run

1. Open `Weesp.Rproj` in RStudio. This sets the working directory to this folder.
2. Open `Weesp_analysis.R` and source it (Ctrl+Shift+S), or run it section by section
   (document outline: Ctrl+Shift+O).

The one script produces every figure, table and p-value of the paper and the supplement,
in the order of the paper. It reads its input from `data/` and writes to `output/`:

| Folder | Contents |
|---|---|
| `output/figures_pdf/` | Figures 1–7 and S1–S5 as vector PDF (Helvetica) |
| `output/figures_png/` | The same figures as PNG, 300 dpi |
| `output/tables/` | Table 1, Tables S1–S15 and the supporting tables, as CSV; `sessionInfo_Weesp.txt` |
| `output/saved_results/` | Each completed surrogate run. A run that is saved is loaded instead of recomputed; delete a file, or set `use_saved_results <- FALSE`, to rerun it |

## Data (`data/`)

| File | Content |
|---|---|
| `AV_diagnostic_data.csv`, `AV_constant40.csv` | The std050, vsc101 and constant40 chronologies with sample depth |
| `stal_scot.csv` | Uamh an Tartair composite growth rate (Baker et al., 2015) |
| `sawtooth2020.txt` | South Sawtooth Lake AMV index (Lapointe et al., 2020), as archived at NOAA |
| `palgeo.jpg` | Map of Figure 1, modified after Vos et al. (2018) |
| `AV_objects_07.rds` (optional) | The detrending variants of S6, built from the raw ring widths by J. Bazelmans' step scripts. Without it, that one section is skipped |

## Surrogate counts and seeds

Set at the top of the script. The published values use 10,000 surrogates for Table 1,
2,000 surrogate fields for the persistence test, 20,000 for the band-pass correlations,
500 for the step test (Figure S5), 1,000 for Table S13, and 2,000 surrogates and 1,000
fields for the supplementary suite at ω₀ = 6 (Tables S1–S4, S6–S8). Seeds: 1337 for
the main analysis; Table S13 (seed 2013), the supplementary suite and the detrending
variants (seed 20092026) each run on their own random number stream, so every section
returns the same values whether or not the others have been run. A complete run at the
published counts takes several hours.

## Packages

`astrochron`, `multitaper` and `readxl` are required; `scico` (colour maps) and `jpeg`
(Figure 1) are optional. Only base graphics and the standard `pdf()` and `png()` devices
are used.
