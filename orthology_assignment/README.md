# Human–Mouse Orthology Assignment

This folder contains the files used to classify human–mouse orthology 
relationships.

## Files

* `Human_to_Mouse_Orthologs_HMD_HumanPhenotype.rpt`
  Input human–mouse orthology table obtained from the Mouse Genome 
Informatics database (Jackson Lab).

* `orthology_assignment.R`
  R script that classifies each ortholog pair as:

  * one-to-one
  * one-to-many (human-to-mouse)
  * one-to-many (mouse-to-human)
  * many-to-many

* `orthologs_with_classification.csv`
  Output table containing the original orthology information and the 
assigned orthology classification.

## Usage

Open `orthology_assignment.R` in R or RStudio and run the script from this 
folder.

The script requires the R package:

```r
install.packages("dplyr")
```

Running the script will generate or overwrite:

```text
orthologs_with_classification.csv
```

