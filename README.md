# Effectiveness of PRS models across populations for cardiovascular disease 🧬


This project aims to investigate the generalizability of PRS models across various populations for cardiovascular diseases. Various models trained on European data are compared against GWAS summary statistics for European, East Asian, African, and American populations. These models and association tests are done using the [`FUSION`](http://gusevlab.org/projects/fusion/) software. If you're interested in learning more about this project, check out our [website](https://y00628.github.io/equitable-prs-cardiovascular/) and [report](https://drive.google.com/file/d/1zu4fIurSg-mX8VTd10HnFDjyd7nwPdN7/view?usp=drive_link)!




## Setup

1. This code requires Python and R to be added to `PATH`.

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

        * Try the following if you encounter an error: 
        `Rscript -e "Rscript -e "devtools::install_github('gabraham/plink2R/plink2R')"`

    * `Rscript install.packages(c('optparse','RColorBrewer'))`

        * Again you may need to try: 
        `Rscript -e "install.packages(c('optparse','RColorBrewer'), repos='http://cran.rstudio.com/')"`

    * `Rscript install.packages(c('glmnet','methods'))`

        * Again you may need to try: 
        `Rscript -e "install.packages(c('glmnet','methods'), repos='http://cran.rstudio.com/')"`

1. Install GCTA (Genome-wide Complex Trait Analysis) if needed
   
 * TWAS FUSION comes with `gcta_nr_robust` (for non-Windows) and `gcta_nr_robust.exe` (for Windows). If it's not working, follow these steps:
   
     * Download the corresponding GCTA version from [here](https://yanglab.westlake.edu.cn/software/gcta/#Download).
       
     * Unzip the downloaded folder.
        
     * Move the executable file into the `src/fusion_twas-master` directory and rename it to `gcta_nr_robust`.
         * e.g., on macOS ARM64, the executable file is named `gcta64`


## Running the code
1. The code can be run using:

    * Linux/macOS: `python3 run.py <arg1> <arg2> ...`

    * Windows: `py run.py <arg1> <arg2> ...`

1. A list of arguments for the program can be listed using the `help` argument:
    
    *Linux/macOS: `python3 run.py help`

    * Windows: `py run.py help`

1. The `data` argument downloads the GBMI GWAS summary statistics for various
populations (European, East Asian, American, African) and processes the
data for use in [FUSION.compute_weights.R](http://gusevlab.org/projects/fusion/#compute-your-own-predictive-models).

    * The GWAS data may take a few minutes to download, depending on the speed of your internet connection, as the files are ~2GB total. These files will be saved in 2 locations: `data/raw/gwas/` and `data/gwas`.

1. The `compute_weights` argument should be used once all data has been 
downloaded and processed using the `data` argument. This will produce model weights required for
our [FUSION.assoc_test.R](http://gusevlab.org/projects/fusion/#typical-analysis-and-output), which will be used to evaluate a European population
based model against all the populations mentioned above. The weights will be stored in `data/weights`

1. The `association_tests` argument should be used once all the model weights have been computed, and will store the test outputs for each population in `data/gene_associations`.

1. The `generate_plots` argument should be used once all other arguments have been run. This will generate plots related to our analysis of the association tests run for each population, which will be stored in the `plots` directory.

