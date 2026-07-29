# Mouse Secretome Annotation

This folder contains the scripts and supporting files used to classify mouse genes as secreted or non-secreted and to add secretome and orthology annotations to the tissue-specific expression mastersheets.

## Input

TPM expression matrix:

```text
gene_expression_matrix_C57_CD1_TPMs.csv
```

## Workflow

Run the scripts in order from the main `Mouse_Transcriptome_SG` directory:

```bash
Rscript Mouse_secretome_annotaton/01_secretome_annotation_dictionary.R
Rscript Mouse_secretome_annotaton/02-annotate_each_gland_secreted_non_secreted_mouse.R
Rscript Mouse_secretome_annotaton/03-adding_orthology.R
```

### 1. Build the secretome dictionary

`01_secretome_annotation_dictionary.R` creates the gene-level secretome annotation dictionary.

Outputs:

```text
Geneid_secretome_dictionary.csv
Geneid_secretome_dictionary.rds
Geneid_secretome_dictionary_not_curated.csv
Geneid_secretome_dictionary_not_curated.rds
```

### 2. Annotate tissue mastersheets

`02-annotate_each_gland_secreted_non_secreted_mouse.R` adds secretome annotations to the PAR, SM, SL, LIV, and PANC expression mastersheets.

The resulting tables include information such as:

* `secreted`
* `Reference`
* `Function`
* `Mean_TPM`
* `Mean_TPM_Male`
* `Mean_TPM_Fem`

### 3. Add orthology annotations

`03-adding_orthology.R` adds mouse–human orthology information to the annotated mastersheets.

## Manual updates

To add or edit a gene function, update:

```text
function_add.csv
```

To manually classify a gene as secreted or non-secreted, update:

```text
mannual_curation.csv
```

After editing either file, rerun the relevant scripts beginning with:

```bash
Rscript Mouse_secretome_annotaton/01_secretome_annotation_dictionary.R
```

## Notes

Male and female expression means are identified from sample names containing `-Mal-` and `-Fem-`.

Known contaminated parotid samples are removed automatically by the annotation scripts.
