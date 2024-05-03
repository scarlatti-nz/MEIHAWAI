import sys
import shutil
import os
output_dir = sys.argv[1]

for dir in os.listdir(output_dir):
    grids = os.listdir(output_dir+dir)
    for grid in grids:
        if 'grid' in grid:
            shutil.rmtree(output_dir+dir+'/'+grid)

shutil.rmtree('logs')