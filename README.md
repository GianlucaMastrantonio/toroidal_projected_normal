# toroidal_projected_normal

This folder contains all the code and data used for the analyses and results for the paper ``*An interpretable family of projected normal distributions and a related copula model for Bayesian analysis of hypertoroidal data*´´


## How to Use

1. **Run Simulations**
   - For Toroidal Projected Normal model:  
     ```r
     source("simulations/tpn_simulations.R")
     ```
   - For Wrapped Cauchy model:  
     ```r
     source("simulations/cwc_simulations.R")
     ```

2. **Analyze Real Data**  
   - Wrapped Cauchy:  
     ```r
     source("realdata/cwc_real.R")
     ```
   - Toroidal Projected normal:  
     ```r
     source("realdata/tpn_real.R")
     ```

3. **Generate Figures and Tables**  
   - To generate correlation plots:  
     ```r
     source("plots and tables/Correlations.R")
     ```
   - Plots for the real data application can be generated after running the two files in the folder **simulations/**:  
     ```r
     source("plots and tables/PlotReal.R")
     ```
   - For simulated data, generate plots after running the two files in the folder **real_data/**:  
     ```r
     source("plots and tables/PlotSim.R")
     ```
   

## Notes

- Set your working directory to the `toroidal_projected_normal` folder before running scripts.
- Results and outputs are saved in the corresponding `output/` directories.
- For any problem and/or questions, email gianluca.mastrantonio@polito.it

## License

See [LICENSE](LICENSE) for details.