from datetime import time
import pandas as pd
import pyogrio as rio
import numpy as np
import os
import sys
import subprocess
import pyarrow as pa
import shutil



output_dir = sys.argv[1]
dirs = os.listdir(output_dir)
for timestep in range(0,26,1):
    agent_tables = []
    catchment_tables = []
    extension_tables = []
    for dir in dirs:
        try:
            agent_csv = pa.ipc.open_file(output_dir+dir+'/agents_'+str(timestep)+'.arrow').read_all()
            agent_tables.append(agent_csv)
            catchment_csv = pa.ipc.open_file(output_dir+dir+'/catchments_'+str(timestep)+'.arrow').read_all()
            catchment_tables.append(catchment_csv)
            if timestep>0:
                extension_csv = pa.ipc.open_file(output_dir+dir+'/extension_'+str(timestep)+'.arrow').read_all()          
                extension_tables.append(extension_csv)
        except:
            dirs.pop(dirs.index(dir))
    combined_agents = pa.concat_tables(agent_tables)
    combined_catchments = pa.concat_tables(catchment_tables)
    
    with pa.OSFile(f'{output_dir}agents_{timestep}.arrow', 'wb') as sink:
        with pa.RecordBatchFileWriter(sink, combined_agents.schema) as writer:
            writer.write_table(combined_agents)
    with pa.OSFile(f'{output_dir}catchments_{timestep}.arrow', 'wb') as sink:
        with pa.RecordBatchFileWriter(sink, combined_catchments.schema) as writer:
            writer.write_table(combined_catchments)
    if timestep > 0:
        combined_extensions = pa.concat_tables(extension_tables)
        with pa.OSFile(f'{output_dir}extension_{timestep}.arrow', 'wb') as sink:
            with pa.RecordBatchFileWriter(sink, combined_extensions.schema) as writer:
                writer.write_table(combined_extensions)
combined = []
for dir in dirs:
    try:
        lu_csv = pa.ipc.open_file(output_dir+dir+'/land_use_init.arrow').read_all()
        combined.append(lu_csv)
    except:
        dirs.pop(dirs.index(dir))
    combined_lu = pa.concat_tables(combined)
    with pa.OSFile(f'{output_dir}land_use_init.arrow', 'wb') as sink:
        with pa.RecordBatchFileWriter(sink, combined_lu.schema) as writer:
            writer.write_table(combined_lu)

