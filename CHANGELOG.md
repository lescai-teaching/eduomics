# eduomics: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.1.0 - Dancing Fibonacci - 2025-11-07

### Credits

Special thanks to the following for their contributions to the release:

- [Francesco Lescai](https://github.com/lescai)
- [Lorenzo Sola](https://github.com/LorenzoS96)
- [Davide Bagordo](https://github.com/DavideBag)
- [Simone Carpanzano](https://github.com/carpanz)
- [Mariangela Santorsola](https://github.com/msantorsola)

### Enhancements & fixes

- [PR #83](https://github.com/lescai-teaching/eduomics/pull/83) - Correct emission in the RNAseq workflow; fix the DNAVALIDATION module; modify and fix the SUBVAR module and subworkflow that calls it; add the simulated gene in the meta of the variant calling workflow; update software dependencies in the RNA-seq workflow
- [PR #91](https://github.com/lescai-teaching/eduomics/pull/91) - Impoved variant and RNA scenario by adding gene argument and modifying the Google Gemini prompt
- [PR #92](https://github.com/lescai-teaching/eduomics/pull/92) - Harmonise the reference format in the SUBVAR module: variant outputs now follow the reference format rather than ClinVar’s. The SUBSET_REFERENCES_TO_TARGET subworkflow was updated accordingly
- [PR #93](https://github.com/lescai-teaching/eduomics/pull/93) - Update the PYCONVERTSIM module to correctly extract the gene name
- [PR #94](https://github.com/lescai-teaching/eduomics/pull/94) - Update the RNA-seq scenario emission from the QUANTIFY_DEANALYSIS_ENRICH_VALIDATE subworkflow
- [PR #95](https://github.com/lescai-teaching/eduomics/pull/95) - Integrate SimuSCoP 2.0.1 across eduomics with reproducible environments

### Software dependencies

| Dependency        | Old version | New version |
| ----------------- | ----------- | ----------- |
| `Biostrings`      | 2.68.1      | 2.74.0      |
| `org.Hs.eg.db`    | 3.17.0      | 3.20.0      |
| `clusterProfiler` | 4.8.1       | 4.14.0      |
| `rtracklayer`     | 1.60.0      | 1.66.0      |
| `annotationDbi`   | 1.62.2      | 1.68.0      |
| `igraph`          | 1.5.1       | 2.1.4       |
| `jsonlite`        | 1.8.7       | 2.0.0       |
| `Matrix`          | 1.6.5       | 1.7.4       |
| `SimuSCoP`        | 1.1.2       | 2.0.1       |

## 1.0.0 - Running Franklin - 2025-07-08

First release of eduomics, created with the [nf-core](https://nf-co.re/) template.

### Credits

Special thanks to the following for their contributions to the release:

- [Francesco Lescai](https://github.com/lescai)
- [Mariangela Santorsola](https://github.com/msantorsola)
- [Lorenzo Sola](https://github.com/LorenzoS96)
- [Davide Bagordo](https://github.com/DavideBag)
- [Simone Carpanzano](https://github.com/carpanz)

## v1.0.0dev - [2025-03-16]

Initial development of eduomics, created with the [nf-core](https://nf-co.re/) template.
