
### SETTING UP YOUR WORKING DIRECTORY

setwd("your own WD")   
getwd()          
list.files()     

source("install_dependencies.R")   # Only the first time

### Running the entire pipeline
source("run_all.R")

### SENSITIVITY ANALYSIS EXCLUDING ITEM 34 with different marginal distribution of distractors

Sys.setenv(MATRIKS_EXCLUDE_ITEMS = "34", MATRIKS_OUTPUT_DIR = "outputs_sensitivity")
source("R/00_setup.R")
source("R/03_matching_and_models.R"); source("R/04_posterior_contrasts.R")
         
