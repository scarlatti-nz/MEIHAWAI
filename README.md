### Repository for the MEIHAWAI (Measuring the Effects of Interventions with Heterogeneous Agents on Water using Agent Interactions) model ###

MEIHAWAI is a national scale, hectare resolution simulation model of producer behaviour. It is designed to simulate the effects of policy interventions on land use, water quality, and farm profitability in Aotearoa New Zealand. The model is based on the principles of agent-based modelling, and simulates the behaviour of individual farmers in response to policy interventions, market conditions, and environmental factors. The model is designed to be used as a tool for policy analysis, and can be used to simulate the effects of a wide range of policy interventions on land use, water quality, and farm profitability. You can read more about the model in the paper linked below.

LINK TO PAPER HERE

This project is distributed under the Creative Commons Attribution-ShareAlike 4.0 International License. See LICENSE.md for details.

### How do I get set up? ###

#### Summary of set up ####
MEIHAWAI is written in [Julia](https://julialang.org) (model engine) and [Python](https://www.python.org) (plotting and data processing utilities), and uses the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system to tie it all together. Since the model is designed to simulate the entirety of Aotearoa New Zealand's productive land with hectare-scale resolution, it makes use of high-level paralellisation, optimised to run on a HPC cluster such as [NeSI](https://www.nesi.org.nz/). The national scale analysis, across the four scenarios specified in the RunEverythingConfig.yaml file, should take around 800 CPU hours to run. Running the workflow on a single machine is possible, but will take several weeks to complete at national scale, so we recommend using a HPC cluster if possibl, or restricting the analysis to a much smaller region of interest.

#### Downloads ####

The model's input data is provided in the form of 783 0.2° x 0.2° geopackage files, available for download from LINK. These files are used to initialise the model and provide the spatial context for the simulation. The model will read all the files in the directory specified in the RunEverythingConfig.yaml file, so to run a smaller scope simulation, simply remove the input files you don't need from the directory.

A single grid cell is provided in the data directory for testing purposes, which can be used to check that the model is running correctly on your system and test alterations to the model.

You will also need to download:

* [Julia](https://julialang.org/downloads/) - the model was developed using Julia 1.9.1, but should work with any version of Julia >=1.9.0

* [Python](https://www.python.org/downloads/) - version 3.8. We recommend using the [Conda](https://docs.conda.io/en/latest/miniconda.html) package manager to install Python and the required packages, but any Python 3.8 installation should work.

Once you've installed Julia and Python, you can install the required packages by running the following commands in the terminal:

```pip install -r requirements.txt```

```julia --project=/./ -e 'using Pkg; Pkg.instantiate()'```


#### Running the model ####

The entire model process, from inputs to report-quality figures, is managed by the snakemake workflow management system. To run the model, simply run the following command in the terminal with your python environment active:

```snakemake -c8 --snakefile RunEverything```

By default, this runs four modelling scenarios in parallel, using 8 cores (the -c8 flag, which can be adjusted to suit your system). The RunEverythingConfig.yaml specifies the modelling parameters used in each scenario, as well as the directory in which model outputs are saved, and a different directory for storing the final report-quality figures. Detailed instructions on how to use the config file are specified in comments in the file itself. 

To run the model on a HPC cluster, you can create a cluster-config.yaml file from the cluster-config.yaml.tmpl provided, and run the following command:

```snakemake --jobs 800 --printshellcmds --rerun-incomplete --snakefile RunEverything --cluster-config cluster-config.yaml --cluster "sbatch --account={cluster.account} --partition={cluster.partition} --mem={cluster.mem} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --time={cluster.time} --hint={cluster.hint} --output={cluster.output} --error={cluster.error}"```

#### Creating your own scenarios ####

You can create your own scenarios by editing the RunEverythingConfig.yaml file. You can do this by changing parameters under the existing run IDs, which will feed through to the plotting utility, or by adding new run IDs and specifying the parameters you want to change. In this case, your new scenario will not be automatically added to the plotting utility, so you will need to add it manually or perform your own analysis using the outputs in the output directory. 

Snakemake will not re-run the model for a specific run ID if the output files already exist, so if you want to re-run the model with different parameters, you will need to delete or move the output files for that run ID.

You can also edit lower level parameters in the inputs spreadsheet data/Model Inputs.xlxs, or even in the code itself. We bear no responsibility for the consequences of doing this, however, and recommend you only do so if you are confident in your understanding of the model and its inputs. 

#### Auxiliary scripts

We include a number of auxiliary scripts used in testing, development, and calibration of the model in the `calibration` and `utils` folders. These scripts are not necessary for running the model, but may be useful for understanding how the model works, or for developing new features. Some of these depend on large datafiles not included in this repository, in these cases a note is included in the script itself describing how to obtain the necessary data. We do not expect these scripts to be useful to most users, but they are included for completeness.

### Who do I talk to? ###

Please direct all feedback to Kenny Bell (kenny.bell@scarlatti.co.nz), or open an issue on the GitHub repository.
