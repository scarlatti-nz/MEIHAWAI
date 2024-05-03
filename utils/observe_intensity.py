import pandas as pd
import numpy as np

baseline = pd.read_feather("../nesi_downloads/agents_0.arrow")
res = []
for lu in pd.unique(baseline['land_use_name']):
    print(lu)
    baseline_results = baseline[baseline['land_use_name'] == lu]
    # for typology in pd.unique(baseline_results['max_nitrogen_loss']):
    #     parcels = baseline_results[(baseline_results['land_use_name'] == lu) & (baseline_results['max_nitrogen_loss'] == typology)]
    #     if len(parcels) > 100:
    #         fifth = np.percentile(baseline_results[(baseline_results['land_use_name'] == lu) & (baseline_results['max_nitrogen_loss'] == typology)]['intensity'],5)
    #         median = np.percentile(baseline_results[(baseline_results['land_use_name'] == lu) & (baseline_results['max_nitrogen_loss'] == typology)]['intensity'],50)
    #         ninetyfifth = np.percentile(baseline_results[(baseline_results['land_use_name'] == lu) & (baseline_results['max_nitrogen_loss'] == typology)]['intensity'],95)
    #     else:
    parcels = baseline_results[(baseline_results['land_use_name'] == lu)]
    typology = "all"
    fifth = np.percentile(baseline_results[(baseline_results['land_use_name'] == lu)]['intensity'],5)
    median = np.percentile(baseline_results[(baseline_results['land_use_name'] == lu)]['intensity'],50)
    ninetyfifth = np.percentile(baseline_results[(baseline_results['land_use_name'] == lu)]['intensity'],95)
    res.append([lu,typology,fifth,median,ninetyfifth])
df = pd.DataFrame(res,columns=['land_use_name','max_nitrogen_loss','fifth','median','ninetyfifth'])
df.to_csv("calibration/observed_intensities.csv")
       