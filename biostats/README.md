# Biostats R environment

See the [portable install and run instructions](../README.md). From **your checkout's** `biostats/` directory, install the packages pinned by `renv.lock` once:

```bash
rig run -r 4.6.1 -f bootstrap.R
```

Create a batch-owned **study** recipe with the Rust TUI, then execute the **validate**, **plan** and **run** commands it prints. The study preprocesses all biological samples and seals their shared PCA/PCoA. To add independent biological comparisons later, reopen the TUI with `--batch <batch> --study <full-study-run_hash>` and execute the same three commands for each new comparison recipe. Comparison analyses reuse the sealed study without changing its ordinations. `run` publishes the dashboard and factual README under `<batch>/results/stats/`; an existing completed hash is republished without recomputing it. Input, code, package environment and output checksums remain recorded with every seal.

The saved PLS-DA validation status means the configured resampling and permutations completed, not that the classifier discriminates reliably. Full-fit loadings and VIP scores are descriptive, not held-out evidence of individual feature importance. Review group sizes, design, permutation results and uncertainty before scientific interpretation.
