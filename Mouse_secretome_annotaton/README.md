**Mouse Secretome & Mastersheets**

TPMs in: gene_expression_matrix_C57_CD1_TPMs.csv

**Build the dictionary**
Rscript secretome_annotation_dictionary.R

**Build per-gland mastersheets (PAR, SM, SL, LIV, PANC)**
Rscript annotate_each_gland_secreted_non_secreted_mouse.R

**How to update Function**
Open mouse_proteome/function_add.csv. Add (or edit) the row for your Geneid with Function (and optional Reference).
Re-run:
Rscript secretome_annotation_dictionary.R
Rscript annotate_each_gland_secreted_non_secreted_mouse.R
Rscript adding_orthology_to_mouse_sheets_update.R

**How to set secreted manually**
Open mouse_proteome/mannual_curation.csv.
Add (or edit) the row for your Geneid with secreted = manual_annotation
(optional: add Reference / Function).
Re-run:
Rscript secretome_annotation_dictionary.R
Rscript annotate_each_gland_secreted_non_secreted_mouse.R

**What you get**
Dictionary: mouse_proteome/Geneid_secretome_dictionary.csv (+ .rds)
Mastersheets: {PAR,SM,SL,LIV,PANC}_mastersheet_TPMs_annotated.csv
(includes secreted, Reference, Function, Mean_TPM, Mean_TPM_Male, Mean_TPM_Fem)

**Notes**
Sex means use column names with -Mal- / -Fem-.
Known contaminated PAR columns are dropped automatically.
