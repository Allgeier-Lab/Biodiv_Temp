## 01_Scripts
Contains all scripts needed to reproduce results and figures. 

### project_code.R
- Script used to generate submission script for High-Performance Computing (HPC).
- Establishes input variables for the HPC submission.
- Declares function for HPC submission that:
  - Runs the experimental individual-based model.
  - Processes the resulting model run, calculating primary production (ecosystem, reef, and open seagrass), average total fish biomass,
 and production per unit fish biomass (ecosystem, reef, and open seagrass).
  - Returns these cumulative metrics for the final timestep of the model run.

### plots.R
  - Uses data generated in project_code.R to create plots.
  - Function titles are informative of the resulting plot (i.e. "plot_production_H1()" will plot primary production values in response to 
the first hypothesis).
