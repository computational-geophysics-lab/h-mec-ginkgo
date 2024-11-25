# Making plots with the julia script:
The julia script requires julia and the packages "Plots" and "DelimitedFiles".
To generate plots from a experiment, open a terminal and enter the julia REPL
```bash
> julia
```
Then do:
```julia
include("GeneratePlots.jl")
create_plots("path/to/your/experiments/h-mec-rsf-XXX-vX/output_data")
```
The script will create a folder "plots" inside the output_data folder that contains all plots corresponding to each of the EVO_xxx.txt files.

For now it only works if the grid size is 401x51.
