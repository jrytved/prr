# prr
Proteomics data processing &amp; analysis in R

## Package outline
This package contains functions that i generally will reuse many times when analyzing proteomics data. Documentation is handled by roxygen2. All code is contained in .R files that are logically grouped. 

## Package structure

### diann_phospho.R
Any functions related specifically to handling phosphoproteomic search outputs from DIA-NN.
### diann.R
Any functions related to handling arbitrary search outputs from DIA-NN.
### peptides.R
Any functions related to peptides and computation of their physicochemical properties.
### plotting.R
Any functions related to plotting data, mainly with ggplot.
### sequence.R
Any functions related to manipulating protein sequences.
### setops.R
Any functions related to set-wise operation operating on for example between sets of proteins.
### tables.R
Any functions related to displaying data as tables, mainly with gt.