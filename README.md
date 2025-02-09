# Effectiveness of PRS models across populations for cardiovascular disease







## Setup

1. This code requires python and R to be added to `PATH`.

1. Clone this repository using the following command: 

    `git clone https://github.com/y00628/equitable-prs-cardiovascular.git`

    If you are having problems due to the size of the repo causing early EOF try these steps in your terminal:
    
    1. `git config --global http.postBuffer 524288000`
    
    2. `git clone --config core.compression=0 https://github.com/y00628/equitable-prs-cardiovascular.git`

1. Navigate to the repositories `src/` folder using:

    `cd path/to/equitable-prs-cardiovascular/src`

1. Install the required python packages by running the following command: 

    `pip install -r requirements.txt` 

1. Install the R packages required by [`TWAS FUSION`](http://gusevlab.org/projects/fusion/) by running the following commands: 

    * `Rscript devtools::install_github("gabraham/plink2R/plink2R")`

    * `Rscript install.packages(c('optparse','RColorBrewer'))`

    * `Rscript install.packages(c('glmnet','methods'))`

1. The code can be run using:
    `py run.py <arg1> <arg2> ...`

1. A list of arguments for the program can be listed using the `help` argument:
    `py run.py help`

1. The `data` argument downloads the GBMI GWAS summary statistics for various
populations (European, East Asian, American, African, All) and processes the
data for use in [FUSION.compute_weights.R](http://gusevlab.org/projects/fusion/#compute-your-own-predictive-models).

1. The `compute_weights` argument should be used once all data has been 
downloaded and processed. This will produce model weights required for
our [FUSION.assoc_test.R](http://gusevlab.org/projects/fusion/#typical-analysis-and-output), which will be used to evaluate a European population
based model against all the populations mentioned above.