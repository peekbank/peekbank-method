# repository for methodological analyses of Peekbank

[Peekbank](http://peekbank.stanford.edu)

note that some intro writing is in here: https://docs.google.com/document/d/1CCp_0hFNctYwmXsj0meBs-XrdssG56ZFGzMFMx92a5M/edit

## Pipeline

### Step 1: Data loading (`0_get_data.Rmd`)

Fetches data from the Peekbank database (version 2026.1), joins AOI timepoints to administrations, trials, and stimuli, filters to vanilla trials and ages ≤ 60 months, and applies item cleanup using `metadata/included_items_2026.1.csv`. Also fuzzy-joins CDI data to administrations.

Most chunks have `eval=F` and are meant to be run once; re-run when the Peekbank database version changes. Produces:
- `cached_intermediates/0_aoi_data_joined.Rds` — raw joined AOI data
- `cached_intermediates/0_d_aoi.rds` — filtered, cleaned AOI data
- `cached_intermediates/0_cdi_subjects.rds` — CDI data joined to administrations

### Step 2: Analysis scripts (`scripts/`)

Scripts compute reliability (ICC, test-retest) and validity (CDI correlation) across parameter grids. Results are written to `cached_intermediates/`. Scripts 1–3 can be run sequentially via `script_runner.R`; scripts 4, 5, and 8 are run individually (they are heavier).

| Scripts | Description |
|---------|-------------|
| `1_acc_*.R` | Accuracy (window-averaged proportion correct): ICC, CDI correlation, test-retest |
| `2_bc_*.R` | Baseline-corrected accuracy: ICC, CDI correlation, test-retest |
| `3_make_rts.R` | Compute reaction times from RLE-encoded AOI data |
| `3_rt_*.R` | RT: ICC, CDI correlation, test-retest |
| `4_acc_trial_*.R` | Trial-exclusion parameter sweep |
| `4_rt_*.R` | RT trial-count parameter sweep |
| `5_acc_*.R`, `5_rt_*.R` | Simulated downsampling of trials |
| `8_acc.R`, `8_rt.R` | Recommended and alternative settings only (used in paper) |

Helper files sourced by scripts:
- `helper/common.R` — core functions: data prep, ICC, CDI correlation, test-retest
- `helper/params.R` — parameter grid definitions
- `helper/rt_helper.R` — RT computation from RLE data

### Step 3: Paper (`paper/`)

`paper/peekbank-method-paper.Rmd` and `paper/supplement.Rmd` read from `cached_intermediates/` and produce figures and statistics. They source `helper/notebook_helper.R` for plotting and modeling functions.

---

## Not currently in use

The following are kept for reference but are not part of the active pipeline:

- `1_accuracies.Rmd` – `7_order_effects.Rmd`: old exploratory notebooks at the repo root (superseded by `scripts/`)
- `peekbank-method-summary.Rmd`: old summary notebook
- `cluster_old/`: old cluster computing scripts from an earlier version of the analysis
- `data_practices/`: PDF codings of data practices across datasets
- `rsconnect/`: deployment config from an old Shiny/Connect app
