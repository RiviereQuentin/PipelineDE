# PipelineDE
Differential expression analysis pipeline in an R Shiny app

## Installation

Pipeline is an R shiny app that requires the updated version of R (R > 4.0.4), BiocManager and remotes to be installed. 

In R, type the following lines:
```
    if(!require("remotes", quietly = TRUE)){  
        install.packages("remotes")
        }
    if(!require("BiocManager", quietly = TRUE)){  
        install.packages("BiocManager")
        }
 ```
  
Then, you can enter:
```
  options(repos = BiocManager::repositories())
  getOption("repos")
  BiocManager::install("RiviereQuentin/PipelineDE",                     
    dependencies = TRUE,                     
    build_vignettes = TRUE,
    force = TRUE)
````
## Perform DE analysis

Type in the R console the following code and follow the instructions on the user interface:

````
library(PipelineDE)
PipelineDE()
````
