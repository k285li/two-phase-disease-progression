This folder contains the R code and a simulated data example for the paper:
“Two-Phase Designs for Biomarker Studies when Disease Processes are under Intermittent Observation.”

Contents
share.func.R: Shared helper functions used by the simulation and data-example scripts.
2phase-simulation.R: Main simulation script. Running this script reproduces the simulation results reported in Table 2 of the paper.
2phase-data-example.R: A worked example illustrating implementation of pseudo-score residual-dependent subsampling and maximum likelihood estimation using simulated data.

How to run
1. Open R or RStudio and set the working directory to this folder.
2. Run the simulation: Create a folder named `results`. Source or run 2phase-simulation.R to reproduce Table 2.
3. Run the simulated-data example: Source or run 2phase-data-example.R.

Notes
The scripts assume that share.func.R is available in the same directory (or sourced with the correct relative path).
Output files (if any) will be written to the working directory unless otherwise specified within the scripts.

Contact
For questions about the code, please check the github page: https://github.com/k285li/two-phase-disease-progression, or contact the corresponding author of the manuscript.
