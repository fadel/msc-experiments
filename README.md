# Introduction
These scripts are used to reproduce the experiments I ran during my Master's
research. The R scripts are the main experiments, while the python scripts in
the `datasets` directory are used to download and prepare the datasets used.

# R scripts
The main R script is `run.R`. It is configured via `config.R`, which contains
details about which techniques and datasets we are running the experiments on.

# Python scripts
Most scripts should run with a simple `python script.py`. If it requires
additional parameters, it will say so and explain which are those. Whenever that
happens, the parameters used in my experiments are detailed in a `README` file
inside the corresponding directory.
