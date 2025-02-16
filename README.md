# Covariate-Adjusted Functional Data Analysis for SHM

2025-02-13

- [Reproduction details and
  instructions](#reproduction-details-and-instructions)
- [Extract bridge data](#extract-bridge-data)
- [Case study 1 - computational
  results](#case-study-1---computational-results)
- [Case study 1 - plot results](#case-study-1---plot-results)
- [Introduction and misc plots](#introduction-and-misc-plots)
- [Simulation study: model training pipe line - computational
  results](#simulation-study-model-training-pipe-line---computational-results)
- [Simulation study: model training pipe line - plot
  results](#simulation-study-model-training-pipe-line---plot-results)
- [Simulation study: monitoring scheme run length - computational
  results part
  1](#simulation-study-monitoring-scheme-run-length---computational-results-part-1)
- [Simulation study: monitoring scheme run length - computational
  results part
  2](#simulation-study-monitoring-scheme-run-length---computational-results-part-2)
- [Simulation study: monitoring scheme run length - computational
  results part
  3](#simulation-study-monitoring-scheme-run-length---computational-results-part-3)
- [Simulations monitoring scheme run length -
  plots](#simulations-monitoring-scheme-run-length---plots)
- [References](#references)

# Reproduction details and instructions

The folders contain the following files, which have been used to
generate the results and plots of the paper ([Wittenberg et al.
2024](#ref-Wittenberg.etal_2024)). For each file an input and output is
stated, which describes the required input data and the generated output
files. The files are written in the R language and replications can be
achived by rendering the qarto “.qmd” files in the R environment.

To run all codes the following R-packages are required:

``` r
# Package names
packages <- c("benchmarkme", "dplyr", "forcats", "funData", "ggplot2", "gratia", 
              "lme4", "lubridate", "mgcv", "patchwork", "plotly", "purrr", 
              "refund", "R.matlab", "segmented", "spc", "tidyr", "tidyselect", 
              "xtable") # "reticulate" and "png" are not required.
                        # R package "reticulate" can be used to run Python code in 
                        # R and save the 3 dimensional plots. The R package "png"
                        # can be used to to include png files.          

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
```

# Extract bridge data

The file *Extract_KW51_bridge_data_Maes.Lombaerts_2021.qmd* downloads
and extracts the KW51 data set which is freely available from
https://zenodo.org/records/3745914, cf. ([Maes and Lombaert
2020](#ref-Maes.Lombaert_2020)), and for information on the bridge
compare ([Maes and Lombaert 2021](#ref-Maes.Lombaert_2021)).

- Output:
  - data/Dataset_bridge_KW51.RDS (Compiled data set for the KW51 bridge)

# Case study 1 - computational results

The file *Case_study_1_computational_results_KW51_Mode06.qmd* replicates
all computations of the KW51 case study in Section 3 of the main paper.

- Input: Dataset_bridge_KW51.RDS (Compiled data set KW51 bridge)
- Output:
  - data/mode06_gam2_fpcaRes.RDS (Basic model eigenfunctions)
  - data/mode06_gam2d.RDS (Basic model refitted)
  - data/mode06_dtaiB2.RDS (Data including eigenfrequencies)
  - data/mode06_gam4d.RDS (Reduced model refitted)
  - data/mode06_dtaiB4.RDS (Data including eigenfrequencies)
  - data/mode06_gam6d.RDS (Additive model refitted)
  - data/mode06_dtaiB6.RDS (Data including eigenfrequencies)
  - data/mode06_gam7_fpcaRes.RDS (Interactions model eigenfunctions)
  - data/mode06_gam7d.RDS (Interactions model refitted)
  - data/mode06_dtaiB7.RDS (Data including eigenfrequencies)
  - Table 1 (R squared results of all four models)

# Case study 1 - plot results

The file *Case_study_1_plot_results_KW51_Mode06.qmd* replicates all KW51
case study related plots in Section 3 of the main paper.

- Input: All saved output results stored in data folder from file
  Case_study_1_computational_results_KW51_Mode06.qmd
- Output:
  - Figure 4 (Eigenfunctions basic model)
  - Figure 5 (Basic model effects)
  - Figure 6 (Reduced model effects)
  - Figure 7a (Additive model intercept effect)
  - Figure 7b (Additive model covariate effects)
  - Figure 8a (Interaction effect model intercept effect)
  - Figure 8b (Interaction effect model covariate effect)
  - Figure 9 (MEWMA control charts in different setups)

# Introduction and misc plots

The file *Introduction_and_misc_plots.qmd* replicates the introduction
profiles and data generation process plots.

- Input: Dataset_bridge_KW51.RDS (Compiled data set KW51 bridge)
- Output:
  - Figure 1 (KW51 data exemplary response, covariate and error
    profiles)
  - Figure A1 (Artificial data exemplary response, covariate and error
    profiles)

# Simulation study: model training pipe line - computational results

The file
*Simulation_study_model_training_pipe_line_computational_results.qmd*
computes 100 runs with J profiles and saves data to data/simulations
folder.

- Input: -
- Output:
  - generated data files (estimated basic model and eigenfunctions) in a
    specified folder

# Simulation study: model training pipe line - plot results

The file *Simulation_study_model_training_pipe_line_plot_results.qmd*
plots initial functions, simulated data run profiles and their mean.

- Input:
  - data/simulations/Simulation_study_N300_NAOBS0_SEED42_ML.RDS
    (Simulated profiles for J=300)
- Output:
  - Figure A2 (Model training pipeline results on simulated profiles)

# Simulation study: monitoring scheme run length - computational results part 1

The file
*Simulations_study_monitoring_scheme_run_length_computational_results.R*
is an R-script prepared for parallel computing the run-lengths of the
MEWMA control charts.

- Input:
  - fpcaResSim_simulation_runlength.RDS (eigenfunctions for artificial
    data)
  - gam2d_simulation_runlength.RDS (Basic model for artificial data)
  - repeat the script for different lambda settings (0.1, 0.3, 1)
- Output:
  - Run-length files for different MEWMA smoothing parameters

# Simulation study: monitoring scheme run length - computational results part 2

The file *SHEWHART_ARL_calculation_daily.R* is an R-script prepared for
parallel computing the run-lengths of the daily Shewhart control chart
based on a fit of a piecewise linear model where hourly data (residuals)
are averaged to 24h.

- Input:
  - gam2d_simulation_runlength.RDS (Load the exact same data as for the
    MEWMA charts from the gam2d_simulation_runlength_XX.RDS files)
- Output:
  - Run-length files for daily Shewhart charts

# Simulation study: monitoring scheme run length - computational results part 3

The file *SHEWHART_ARL_calculation_hourly.R* is an R-script prepared for
parallel computing the run-lengths of the hourly Shewhart control chart
based on a fit of a piecewise linear model.

- Input:
  - gam2d_simulation_runlength.RDS (Load the exact same data as for the
    MEWMA charts from the gam2d_simulation_runlength_XX.RDS files)
- Output:
  - Run-length files for hourly Shewhart charts

# Simulations monitoring scheme run length - plots

The file *Simulations_monitoring_scheme_run_length_plots.qmd* plots ARL
profiles for different MEWMA smooting parameters.

- Input:
  - data/simulations/ (Average Run Length profiles for lambda=0.1,0.3,1)
- Output:
  - Figure A3 (MEWMA charts Average Run Length performance comparison)

# References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Maes.Lombaert_2020" class="csl-entry">

Maes, K., and G Lombaert. 2020. “<span class="nocase">Monitoring Railway
Bridge KW51 Before, During, and After Retrofitting. V1.0</span>.”
<https://doi.org/10.5281/zenodo.3745914>.

</div>

<div id="ref-Maes.Lombaert_2021" class="csl-entry">

Maes, K., and G. Lombaert. 2021. “Monitoring Railway Bridge KW51 Before,
During, and After Retrofitting.” *Journal of Bridge Engineering* 26 (3):
04721001. <https://doi.org/10.1061/(ASCE)BE.1943-5592.0001668>.

</div>

<div id="ref-Wittenberg.etal_2024" class="csl-entry">

Wittenberg, Philipp, Lizzie Neumann, Alexander Mendler, and Jan
Gertheiss. 2024. “Covariate-Adjusted Functional Data Analysis for
Structural Health Monitoring.” <https://arxiv.org/abs/2408.02106>.

</div>

</div>
