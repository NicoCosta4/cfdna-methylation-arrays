# Reference lists

| File | Content | Source |
|---|---|---|
| `crosshyb_probes.txt` | Cross-reactive probes (450K and EPIC v1): probe sequences that align to non-target regions of the genome (≥ 47 bp match) | ExperimentHub `EH3129`, as used by [DMRcate](https://bioconductor.org/packages/DMRcate/). Compiled from Chen *et al.* (2013) and Pidsley *et al.* (2016) |

Regenerate with:

```bash
docker run --rm -v "$PWD":/project cfdna-methylation-arrays Rscript scripts/fetch_reference_lists.R
```

**References**

- Chen YA, *et al.* Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray. *Epigenetics* 8, 203–209 (2013).
- Pidsley R, *et al.* Critical evaluation of the Illumina MethylationEPIC BeadChip microarray for whole-genome DNA methylation profiling. *Genome Biology* 17, 208 (2016).
