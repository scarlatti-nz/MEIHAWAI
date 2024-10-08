import os
import json
import datetime
import pandas as pd
import random 

configfile: "config.yaml"

run_id = config.get('timestamp',datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")) + 'e' + str(config['extensionScale']) + 'n' + str(config['NPrice']) + 'p' + str(config['PPrice']) + 'sed' + str(config['SedimentPrice']) + 'n2o' + str(config['N2OPrice']) + 'ch4' + str(config['CH4Price']) + 'co2' + str(config['CO2Price']) + 'b' + str(config['enableBolus']) + 'noluc' + str(config['disableLUC']) + 'path'+ str(config['carbonGrowthPath'])

output_dir = config['outputDir']
data_dir = config['dataDir']
all_cells = os.listdir(data_dir)

census = pd.read_csv('utils/summary.csv')
grid_cells = census[(census['grid'].isin(all_cells))]['grid'].tolist()
n_cells = min(len(grid_cells),config.get('nCells',100))
grid_cells = grid_cells[:n_cells]
verbosePlots = "--verbose" if config.get('verbosePlots',False) else ""


print("Running on " + str(len(grid_cells)) + " cells")
rule cleanup:
    input:
        output_dir + run_id + "/plots/initial_land_use.tif"
    shell:
        "python utils/cleanup.py {output_dir} {run_id}"

rule plots:
    input:
        output_dir + run_id + "/agents_25.arrow",
    output:
        plot1 = output_dir + run_id + "/plots/initial_land_use.tif"
    shell:  
        "python utils/generate_plots.py --run {output_dir}{run_id}/ {verbosePlots}  && echo 'plots done!'"


rule join:
    input:
        expand(
            output_dir + run_id + "/{grid_cell}/complete.txt",
            grid_cell = grid_cells
        ),
    output:
        output_dir + run_id + "/agents_25.arrow",
    shell:
        "python utils/aggregate_data.py {output_dir}{run_id}/"


rule main:
    input:
        data_dir + "{grid_cell}/hectares_with_data_final.gpkg",
    output:
        output_folder = output_dir + run_id + "/{grid_cell}/complete.txt",
    shell:
        '''
        julia --project=\".\" main.jl \
        --max-parcels {config[maxParcelsPerCell]} \
        --geodata {data_dir}{wildcards.grid_cell}/hectares_with_data_final.gpkg \
        --output-dir {output_dir}{run_id}/{wildcards.grid_cell}/ \
        --log-file {output_dir}{run_id}/{wildcards.grid_cell}/error.log \
        --extension-scale {config[extensionScale]} \
        --NPrice {config[NPrice]} --PPrice {config[PPrice]} \
        --SedimentPrice {config[SedimentPrice]} \
        --N2OPrice {config[N2OPrice]} \
        --CH4Price {config[CH4Price]} \
        --CO2Price {config[CO2Price]} \
        --sediment-model {config[sedimentModel]} \
        --enable-bolus {config[enableBolus]} \
        --disable-land-use-change {config[disableLUC]} \
        --carbon-growth-path {config[carbonGrowthPath]} \
        --initialise-only {config[initialiseOnly]}
        '''



