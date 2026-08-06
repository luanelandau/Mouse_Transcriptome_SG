# SM sex-DEG genomic clustering analysis

This analysis tests whether sex-associated differentially expressed genes (DEGs) in the mouse submandibular gland occur in 100-kb genomic clusters more often than expected for random sets of genes.

## Inputs

- `generate_bed_SM_sex.R`: joins DESeq2 results to gene coordinates.
- `SM_sex_DEGs.bed`: the 1,615 significant SM sex DEGs (`padj < 0.05` and `|log2FoldChange| > 1`).
- `DEGs.bed`: all genes tested by DESeq2 with available coordinates; this is the permutation background.
- BED columns are chromosome, start, end, and gene ID.

Generate the BED inputs with:

```bash
Rscript generate_bed_SM_sex.R
```

## Observed 100-kb clusters

Run from this directory:

```bash
bedtools sort -i SM_sex_DEGs.bed > SM_sex_DEGs_sorted.bed

bedtools cluster \
  -i SM_sex_DEGs_sorted.bed \
  -d 100000 \
  > SM_sex_DEGs_clusters_100k.bed
```

`bedtools cluster` assigns genes to the same cluster when consecutive intervals are no more than 100 kb apart. Clustering is transitive, so the total span of a cluster can exceed 100 kb.

Summarize the observed number of clusters of each size:

```bash
cut -f5 SM_sex_DEGs_clusters_100k.bed \
  | sort -n \
  | uniq -c \
  | awk '{print $1}' \
  | sort -n \
  | uniq -c \
  | awk 'BEGIN{OFS="\t"} {print $2, $1}' \
  > SM_DEGs_observed_size_distribution_100kb.tsv
```

This headerless file contains `size` and `n_clusters`.

## Permutations

Each of 1,000 permutations samples 1,615 genes without replacement from `DEGs.bed`, then applies the same sorting and 100-kb clustering used for the observed DEGs.

```bash
ALL_GENES="DEGs.bed"
N=1615
N_PERM=1000
DIST=100000
OUTDIR="tmp_perms_100kb"

mkdir -p "$OUTDIR"

for i in $(seq 1 "$N_PERM"); do
  printf -v TAG "%04d" "$i"
  SAMPLE="$OUTDIR/sample_${TAG}.bed"
  CLUST="$OUTDIR/clustered_${TAG}.bed"

  shuf -n "$N" "$ALL_GENES" > "$SAMPLE"
  bedtools sort -i "$SAMPLE" \
    | bedtools cluster -d "$DIST" \
    > "$CLUST"

  if (( i % 50 == 0 )); then
    echo "Finished permutation $i"
  fi
done
```

Summarize all permutation clusters:

```bash
printf 'perm\tsize\tn_clusters\n' > size_distribution_100kb.tsv

for f in tmp_perms_100kb/clustered_*.bed; do
  perm=$(basename "$f" .bed)
  perm=${perm#clustered_}

  cut -f5 "$f" \
    | sort -n \
    | uniq -c \
    | awk '{print $1}' \
    | sort -n \
    | uniq -c \
    | awk -v p="$perm" 'BEGIN{OFS="\t"} {print p, $2, $1}' \
    >> size_distribution_100kb.tsv
done
```

The output columns are `perm`, `size`, and `n_clusters`.

## Statistical analysis

Run:

```bash
Rscript Permutations_analysis.R
```

For each exact cluster size, the script compares the observed number of clusters with the distribution across the 1,000 randomized gene sets. The one-sided empirical p-value is:

```text
(number of permutations with n_clusters >= observed + 1) / (1000 + 1)
```

The script also applies false-discovery-rate (FDR) correction across cluster sizes. Missing cluster sizes in a permutation are treated as zero.

## Main outputs

Outputs are written to `100kb_analysis/`:

- `cluster_size_empirical_pvalues_100kb.csv`: observed counts, null mean and SD, empirical p-values, and FDR-adjusted p-values.
- `cluster_size_permutation_cutoffs_100kb.csv`: 95th- and 99th-percentile null cutoffs by cluster size.
- `cluster_size_key_tests_100kb.txt`: tests for clusters of exactly six genes and at least six genes.
- `cluster_size_distribution_100kb.png/.svg`: permutation distributions with observed values.
- `cluster_size_expected_vs_observed_100kb.png/.svg`: observed and expected lines with an expected +/-1 SD band.

The older `SM_sex_DEGs_clusters_1mb_update.bed` and `tmp_perms/` files are from the earlier 1-Mb analysis and are not inputs to the current 100-kb analysis.
