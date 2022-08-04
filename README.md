# Project

## Title

Reduced basal ganglia N-acetylaspartate linked to working memory deficits in late-onset versus average-onset first episode schizophrenia-spectrum disorder

## Folder Structure

```bash
├── mrs_sz_ageofonset.Rproj
├── README.md
├── data
│   ├── processed
│   └── raw
├── outputs
│   ├── figs
│   └── tables
├── scripts
│   ├── 01_descriptive.R
│   ├── 02_inferences.R
└── src
    └── R
        └── stat.R
```

## Scripts

```base
# run and get descriptives of the samples
Rscript scripts/01_descriptive.R
# run and anwser the research questions
Rscript scripts/02_inferences.R

# src folder contains source codes for functions and helpers
```
