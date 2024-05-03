import os
import time

def runner_func(co2price):
    os.system("snakemake -c6 --config CO2Price=" + str(co2price))
    dir = "outputs/2024_03_22_18_00_00e0.0n0.0p0.0sed0.0n2o0.0ch40.0co2"+str(co2price)+"bfalsenolucfalsepathparabolic/plots/"
    
    try:
        for i in range(26):
            os.remove(dir[:-6] + "agents_" + str(i) + ".arrow")
            os.remove(dir[:-6] + "catchments_" + str(i) + ".arrow")
            if i>0:
                os.remove(dir[:-6] + "extension_" + str(i) + ".arrow")
    except:
        pass

time.sleep(5000)
for co2price in range(190,205,5):
    runner_func(co2price)