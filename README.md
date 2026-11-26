# Temporal coordination of joint calling in water rails

Code and results of the broad-scale vocal coordination analysis of water rail (*Rallus aquaticus*) pairs.

**Website:** <https://marce10.github.io/waterrail_vocal_coordination/>

**Authors:** Karolina Kuczkowska, David M. Logue, [Marcelo Araya-Salas](https://marce10.github.io/), Jan Jedlikowski

## Overview

We tested whether the call series of water rail pair members occur closer in time and overlap more often than expected by chance. Callers could not be identified, so overlapping series were assigned to different birds and non-overlapping series were assigned at random under nine pair-contribution scenarios. Each assignment was compared with randomizations of the call sequence that keep each bird's calling rate.

## Repository structure

| Path | Contents |
|:-----|:-----|
| `index.qmd`, `_quarto.yml`, `styles.css`, `_figure-download.html` | Website pages and configuration |
| `scripts/coordination_analysis.qmd` | Assignment of call series, randomization tests and coordination scores |
| `scripts/null_simulation.qmd` | Validation of the procedure with simulated pairs |
| `scripts/test_coordination_impute.R` | Modified version of `warbleR::test_coordination()` that handles unknown caller identity |
| `output/` | Figures |
| `_freeze/` | Stored results of the report, used to build the website |
| `.github/workflows/publish.yml` | Website deployment to GitHub Pages |

The raw annotations (`data/`) and the manuscript are not in this repository.

## Updating the website

Pushing to `master` publishes the website. The workflow only uses the stored results in `_freeze/`, so it needs neither R nor the data. After changing the report, run `quarto render` in the project root and commit `_freeze/` and `output/`. This needs the local `data/` folder.
