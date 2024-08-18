import seaborn as sb
import pandas as pd
import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scipy.stats as stats
import random
import os
import argparse
import numpy as np
from seaborn import categorical
from sklearn.metrics import confusion_matrix 
import rasterio
from rasterio.features import rasterize
import geopandas as gpd
from shapely.geometry import Point

## initialise
parser = argparse.ArgumentParser()
parser.add_argument("--run",help="run id of the model run to plot",type=str)
parser.add_argument("--verbose",help="increase output verbosity",action="store_true")
args = parser.parse_args()

output_dir = parser.parse_args().run

try:
    os.mkdir(f'{output_dir}plots')
except:
    print("plots directory already exists")


sb.set_palette("bright")
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 13


### run
def import_data():
    

    full_df  = pd.DataFrame()
    catchments = pd.DataFrame()
    extension = pd.DataFrame()
    for i in range(0,26,1):
        df = pd.read_feather(f"{output_dir}agents_{i}.arrow")
        df["Timestep"] = i
        full_df = pd.concat([full_df,df])
        
        catchment_df = pd.read_feather(f"{output_dir}catchments_{i}.arrow")
        catchment_df["Timestep"] = i
        catchments = pd.concat([catchments,catchment_df])
        
        if i > 0:
            extension_df = pd.read_feather(f"{output_dir}extension_{i}.arrow")      
            extension_df["Timestep"] = i
            extension = pd.concat([extension,extension_df])
        
    return full_df, catchments, extension

data, catchment_data, extension = import_data()

lu_to_cat = {"apple":"Horticulture/perennial",
            "avocado":"Horticulture/perennial",
            "blueberry-covered":"Horticulture/perennial",
            "cherry":"Horticulture/perennial",
            "kiwifruit-gold":"Horticulture/perennial",
            "kiwifruit-green":"Horticulture/perennial",
            "maize-grain":"Arable/annual",
            "onions":"Arable/annual",
            "peas-vining":"Arable/annual",
            "potatoes":"Arable/annual",
            "wheat":"Arable/annual",
            "wine-grape-pinot-noir":"Horticulture/perennial",
            "wine-grape-sauvignon-blanc":"Horticulture/perennial",
            "dairy":"Dairy",
            "dairy-irrigated":"Dairy",
            "sheep-and-beef":"Sheep and beef",
            "dairy-bolus":"Dairy",
            "sheep-and-beef-bolus":"Sheep and beef",
            "production-forestry":"Forestry",
            "carbon-forestry":"Forestry",
            "native-forestry":"Forestry",
            "fallow":"Fallow/bare land"}


catchment_data['Primary land use category'] = catchment_data['primary_land_use'].apply(lambda x: lu_to_cat[x])

if args.verbose:
    catchment_data.to_csv(f'{output_dir}plots/catchment_data.csv')

extension = extension.groupby(["Timestep","funding","run_this_timestep"]).sum().reset_index()
extension.to_csv(f'{output_dir}plots/extension.csv')
# retirements = []
# for Timestep in range(1,20):
#     retired = 0
#     for agent in data[(data["Timestep"]==Timestep)]["uid"].unique():
#         age = data[(data["uid"]==agent) & (data["Timestep"]==Timestep)]["age"].values[0]

#         if data[(data["Timestep"]==Timestep+1)&(data["uid"]==agent)]["age"].values[0] != age+1:
#             retired+=1
#     retirements.append([Timestep,retired])
# df = pd.DataFrame(retirements,columns=["Timestep","Number of retirements"])
# sb.barplot(data=df,x="Timestep",y="Number of retirements")
# plt.savefig(f'{output_dir}plots/retirements.png')


# data = data[data["land_use_name"] != "cherry"]
# data = data[data["land_use_name"] != "blueberry-covered"]
# data = data[data["land_use_name"] != "kiwifruit-gold"]





data.sort_values(by=["Timestep","land_use_name"],inplace=True)
data['production_capability'] = data['plant_and_animal_capability']**0.5 * data['people_capability']**0.25 * data['business_capability']**0.25
#rename Environmental sensitivity to Environmental sensitvity (h)
data = data.rename(columns={
                            # "timestep":"Timestep",
                            "area_ha":"Area (ha)",
                            "intensity":"Intensity",
                            "nl_mitigation":"Nutrient Loss Mitigation",
                            "max_revenue":"Max Revenue ($/ha/yr)",
                            # "profitability":"Profit ANPV ($/ha/yr)",
                            "profit_npv_per_ha":"Profit ANPV ($/ha/yr)",
                            "total_yearly_profit":"Yearly profit ($)",
                            "nitrogen_loss":"Nitrogen loss (kg/ha/yr)",
                            "phosphorous_loss":"Phosphorous loss (kg/ha/yr)",
                            "sediment_loss":"Sediment loss (t/ha/yr)",
                            "methane_emissions":"Methane Emissions (kg CO2e/ha/yr)",
                            "nitrous_oxide_emissions":"Nitrous Oxide Emissions (kg CO2e/ha/yr)",
                            "utility":"Utility ANPV ($/ha equivalent)",
                            "land_use_name":"Land use",
                            "years_in_current_land_use":"Time in current land use (yr)",
                            "nl_mitigation_capability":"Nutrient Loss Mitigation Capability",
                            "production_capability":"Production Capability",
                            "deviation_parameter":"Production Capability Error",
                            "time_spent_on_external_resources":"Learning hours delivered",
                            },errors="raise")
data["GHG Emissions (kg CO2e/ha/yr)"] = data["Methane Emissions (kg CO2e/ha/yr)"] + data["Nitrous Oxide Emissions (kg CO2e/ha/yr)"]


data["Total Methane Emissions (kg CO2e/yr)"] = data["Methane Emissions (kg CO2e/ha/yr)"] * data["Area (ha)"]
data["Total Nitrous Oxide Emissions (kg CO2e/yr)"] = data["Nitrous Oxide Emissions (kg CO2e/ha/yr)"] * data["Area (ha)"]
data["Total Nitrogen loss (kg/yr)"] = data["Nitrogen loss (kg/ha/yr)"] * data["Area (ha)"]
data["Total Sediment loss (t/yr)"] = data["Sediment loss (t/ha/yr)"] * data["Area (ha)"]
data["Total Phosphorous loss (kg/yr)"] = data["Phosphorous loss (kg/ha/yr)"] * data["Area (ha)"]

data["subsidy spent"] = data.apply(lambda x: x['Area (ha)'] *  min(x['nl_mitigation_subsidy'], x["Nutrient Loss Mitigation"] * (200 if x["Land use"] == "sheep-and-beef" else 1000)),axis=1)
if args.verbose:
    data.to_csv(f'{output_dir}plots/full_data.csv')


data["Land use category"] = data["Land use"].apply(lambda x: lu_to_cat[x])

sums = data.groupby(["Timestep","Land use"]).sum().reset_index()
averages = data.groupby(["Timestep","Land use"]).mean().reset_index()
averages.to_csv(f'{output_dir}plots/outcomes_average.csv')
sums.to_csv(f'{output_dir}plots/outcomes_aggregated.csv')

sb.relplot(data=sums,x="Timestep",y="Yearly profit ($)",hue="Land use",kind="line")
plt.savefig(f'{output_dir}plots/profitability.png')
plt.close()

aggregates = data.groupby(["Timestep","Land use category"]).sum().reset_index()

for outcome in ["Area (ha)","Yearly profit ($)","Total Methane Emissions (kg CO2e/yr)","Total Nitrous Oxide Emissions (kg CO2e/yr)"]:
    # aggregates.set_index("Timestep").plot(y=outcome,kind="bar",stacked=True)
    sb.relplot(data=aggregates,x="Timestep",y=outcome,hue="Land use category",kind="line")
    plt.savefig(f'{output_dir}plots/{outcome.replace("/"," per ")} by land use.png')
    plt.close()

even_more_aggregated = aggregates.groupby(["Timestep"]).sum().reset_index()
for outcome in ["Yearly profit ($)","Total Methane Emissions (kg CO2e/yr)","Total Nitrous Oxide Emissions (kg CO2e/yr)"]:
    sb.relplot(data=even_more_aggregated,x="Timestep",y=outcome,kind="line")
    plt.savefig(f'{output_dir}plots/{outcome.replace("/"," per ")}.png')
    plt.close()



for lu in pd.unique(data["Land use"]):
    filtered = data[(data["Land use"]==lu)]
    plot = sb.displot(data=filtered[(filtered["Timestep"]==0)],x="Profit ANPV ($/ha/yr)",bins=25)
    plt.savefig(f'{output_dir}plots/profit_hist_{lu}.png')
    plt.close()
    plot = sb.displot(data=filtered[(filtered["Timestep"]==0)],x="Nitrogen loss (kg/ha/yr)",bins=25)
    plt.savefig(f'{output_dir}plots/nitrogen_loss_hist_{lu}.png')
    plt.close()
    plot = sb.displot(data=filtered[(filtered["Timestep"]==0)],x="Phosphorous loss (kg/ha/yr)",bins=25)
    plt.savefig(f'{output_dir}plots/phosphorous_loss_hist_{lu}.png')
    plt.close()
    
    relplot = sb.relplot(data=filtered[(filtered["Timestep"] == 25)],x="Production Capability",y="Profit ANPV ($/ha/yr)",kind="scatter")
    plt.savefig(f'{output_dir}plots/profit_vs_prod_{lu}.png')
    plt.close()
    relplot = sb.relplot(data=filtered[(filtered["Timestep"] == 25)],x="Production Capability Error",y="Profit ANPV ($/ha/yr)",kind="scatter")
    plt.savefig(f'{output_dir}plots/profit_vs_prod_error_{lu}.png')
    plt.close()

unique_agents = data.drop_duplicates(subset=["uid","Land use","Timestep"]).reset_index()
cap_dist_init = pd.DataFrame(columns = ['Land use','mean','stddev'])
cap_dist_init['Land use'] = pd.unique(unique_agents["Land use"])
for lu in pd.unique(unique_agents["Land use"]):
    cap_dist_init.loc[cap_dist_init["Land use"]==lu,'mean'] = unique_agents[(unique_agents["Land use"]==lu) & (unique_agents["Timestep"]==0)]["Production Capability"].mean()
    cap_dist_init.loc[cap_dist_init["Land use"]==lu,'stddev'] = unique_agents[(unique_agents["Land use"]==lu) & (unique_agents["Timestep"]==0)]["Production Capability"].std()
cap_dist_init.loc[len(cap_dist_init.index)] = ["all land uses",unique_agents[(unique_agents["Timestep"]==0)]["Production Capability"].mean(),unique_agents[(unique_agents["Timestep"]==1)]["Production Capability"].std()]
cap_dist_init.to_csv(f'{output_dir}plots/cap_dist_init.csv')


unique_agents_summed = unique_agents.groupby(["Timestep"]).sum().reset_index()
unique_agents_mean = unique_agents.groupby(["Timestep","Land use category"]).mean().reset_index()

sb.lineplot(data=unique_agents_summed,x="Timestep",y="Learning hours delivered")
plt.savefig(f'{output_dir}plots/learning_hours.png')
plt.close()

sb.lineplot(data=unique_agents_mean,x="Timestep",y="Production Capability",hue="Land use category")
plt.savefig(f'{output_dir}plots/mean_prod_cap.png')
plt.close()

sb.lineplot(data=unique_agents_mean,x="Timestep",y="Nutrient Loss Mitigation Capability",hue="Land use category")
plt.savefig(f'{output_dir}plots/mean_nl_cap.png')
plt.close()

#     relplot = sb.regplot(data=filtered[(filtered["Timestep"] == 25) & (filtered["Time in current land use (yr)"] > 1)],x="Production Capability",y="Intensity",lowess=True)
#     plt.savefig(f'{output_dir}plots/prod_vs_intensity_{lu}.png')
#     plt.close()
#     plt.subplots(figsize=(10,6))
#     profit_over_time = sb.relplot(data=filtered[filtered["Timestep"]<=25],x="Time in current land use (yr)",y="Profit ANPV ($/ha/yr)",kind="line")
#     plt.savefig(f'{output_dir}plots/profit_over_time_{lu}.png')
#     plt.close()
#     utility_over_time = sb.relplot(data=filtered[filtered["Timestep"]<=25],x="Time in current land use (yr)",y="Utility ANPV ($/ha equivalent)",kind="line")
#     plt.savefig(f'{output_dir}plots/utility_over_time_{lu}.png')
#     plt.close()



catchment_data['Out of compliance - N'] = catchment_data.apply(lambda x: (x['Load_N'] > x['MAL_N']) & (x['MAL_N'] > 0),axis=1)
catchment_data['Out of compliance - P'] = catchment_data.apply(lambda x: (x['Load_P'] > x['MAL_P']) & (x['MAL_P'] > 0),axis=1)

exempt_segs = pd.unique(catchment_data[catchment_data['has_Sediment_exemption'] == True]['nzseg'])
catchment_data['Out of compliance - Sediment'] = catchment_data.apply(lambda x: (x['Load_Sediment'] > x['MAL_Sediment']) & (x['MAL_Sediment'] > 0) & (x['nzseg'] not in exempt_segs) ,axis=1)  


outcome_time_series = pd.DataFrame(columns = ['Timestep','Out of compliance - N','Out of compliance - P','Out of compliance - Sediment','Methane Emissions','Nitrous Oxide Emissions'])
outcome_time_series['Timestep'] = range(0,26)
for i in range(0,26):
    outcome_time_series.loc[outcome_time_series.Timestep==i,'Out of compliance - N'] = max(1,len(catchment_data[(catchment_data['Timestep']==i) & (catchment_data['Out of compliance - N']==True)]))
    outcome_time_series.loc[outcome_time_series.Timestep==i,'Out of compliance - P'] = max(1,len(catchment_data[(catchment_data['Timestep']==i) & (catchment_data['Out of compliance - P']==True)]))
    outcome_time_series.loc[outcome_time_series.Timestep==i,'Out of compliance - Sediment'] = max(1,len(catchment_data[(catchment_data['Timestep']==i) & (catchment_data['Out of compliance - Sediment']==True)]))


outcome_time_series['Out of compliance - N'] = outcome_time_series['Out of compliance - N'] / outcome_time_series.loc[outcome_time_series["Timestep"]==0 ,'Out of compliance - N'].iloc[0]
outcome_time_series['Out of compliance - P'] = outcome_time_series['Out of compliance - P'] / outcome_time_series.loc[outcome_time_series["Timestep"]==0 ,'Out of compliance - P'].iloc[0]
outcome_time_series['Out of compliance - Sediment'] = outcome_time_series['Out of compliance - Sediment'] / outcome_time_series.loc[outcome_time_series["Timestep"]==0 ,'Out of compliance - Sediment'].iloc[0]



outcome_time_series['Methane Emissions'] = (even_more_aggregated['Total Methane Emissions (kg CO2e/yr)'] / even_more_aggregated.loc[even_more_aggregated["Timestep"]==0,'Total Methane Emissions (kg CO2e/yr)'].iloc[0]).tolist()
outcome_time_series['Nitrous Oxide Emissions'] = (even_more_aggregated['Total Nitrous Oxide Emissions (kg CO2e/yr)'] / even_more_aggregated.loc[even_more_aggregated["Timestep"]==0,'Total Nitrous Oxide Emissions (kg CO2e/yr)'].iloc[0]).tolist()

melted = pd.melt(outcome_time_series,id_vars=["Timestep"])
sb.lineplot(data=melted,x="Timestep",y="value",hue="variable")
plt.savefig(f'{output_dir}plots/outcome_time_series.png')
plt.close()
outcome_time_series.to_csv(f'{output_dir}plots/outcome_time_series.csv')




first_and_last = data[(data["Timestep"]==0) | (data["Timestep"]==25)]

unique_N_catchments = catchment_data[(catchment_data["Timestep"]==0)&(catchment_data['excess_N_percent']>0)]['nzseg'].unique()
unique_P_catchments = catchment_data[(catchment_data["Timestep"]==0)&(catchment_data['excess_P_percent']>0)]['nzseg'].unique()


first_and_last['N_catchment'] = first_and_last["nzseg"].isin(unique_N_catchments)
first_and_last['P_catchment'] = first_and_last["nzseg"].isin(unique_P_catchments)
first_and_last['stayed_in_same_lu'] = first_and_last.loc[first_and_last["Timestep"]==25,"Land use"].eq(first_and_last.loc[first_and_last["Timestep"]==0,"Land use"])
first_and_last.groupby(["Timestep","stayed_in_same_lu","N_catchment","P_catchment"]).sum().reset_index().to_csv(f'{output_dir}plots/changed_lu.csv')


cat_to_class = {"Forestry/other": -1,
                "Forestry":-1,
                "Fallow/bare land":-1, 
                "Dairy": 1,
                "Sheep and beef": 2,
                "Horticulture/perennial": 3,
                "Arable/annual": 4}

lucas_cat_to_class = {"other": -1,
                        "dairy":1,
                        "sheep-and-beef":2,
                        "horticulture":3,
                        "arable":4}



def rasterise_initial_land_use(px):
    df = pd.read_feather(f"{output_dir}land_use_init.arrow")
    df['lucat'] = df["land_use_name"].apply(lambda x: lu_to_cat[x])
    df["class"] = df["lucat"].apply(lambda x: cat_to_class[x])
    df["geometry"] = [Point(xy) for xy in zip(df.x, df.y)]
    gdf = gpd.GeoDataFrame(df, geometry='geometry', crs='EPSG:2193')
    # print(gdf)
    # gdf = gdf.to_crs('EPSG:2193')
    with rasterio.open(f'data/nz_mask.tif') as r:
            nz_mask = r.read(1)

    bbox = [1090035.0758, 4748468.5105, 2089384.2642, 6193899.0394] # bbox in NZTM EPSG:2193

    x_res = (bbox[2] - bbox[0]) / px
    y_res = (bbox[3] - bbox[1]) / px

    # Create an affine transform for the raster
    transform = rasterio.transform.from_origin(bbox[0], bbox[3], x_res, y_res)
    rasterized = rasterize(
        [(geom, value) for geom, value in zip(gdf.geometry, df['class'])],
        out_shape=(px,px),
        fill=0,  # Use a NoData value that does not conflict with your class values
        transform=transform,
        all_touched=True
    )

    mask = (nz_mask == -2)
    rasterized[mask] = -2
    metadata = {
        'driver': 'GTiff',  # GeoTIFF format
        'dtype': 'int32',  # Data type of the raster; adjust as needed
        'nodata': None,  # Define the NoData value; adjust as needed
        'width': px,  # Number of columns of the raster
        'height': px,  # Number of rows of the raster
        'count': 1,  # Number of bands; typically 1 for a classified raster
        'crs': 'EPSG:2193',  # Coordinate Reference System
        'transform': transform
    }

    with rasterio.open(f'{output_dir}/plots/initial_land_use.tif', 'w', **metadata) as dst:
        dst.write(rasterized, 1)
        print("finished write")
    print('exited with block')
    return rasterized
rasterise_initial_land_use(3000)
print('post function call')

# def get_wq_tax_stats():
#     df1 = data[data["Timestep"]==1]
#     nonzero_prices = df1.loc[df1['nitrogen_price']>0,'nitrogen_price']

#     len(nonzero_prices)/len(df1)
#     nonzero_either = df1.loc[(df1['nitrogen_price']>0) | (df1['phosphorous_price']>0)]

#     len(nonzero_either)/len(df1)

#     sum(nonzero_prices)/len(nonzero_prices)

#     nonzero_p_prices = df1.loc[df1['phosphorous_price']>0,'phosphorous_price']
#     sum(nonzero_p_prices)/len(nonzero_p_prices)

#     nonzero_n = df1.loc[df1['nitrogen_price']>0]
#     len(nonzero_n)

#     grouped_to_agents = df1.groupby(['uid']).sum().reset_index()


#     agents_nonzero = grouped_to_agents[grouped_to_agents['wq_pollutant_costs']>0]

#     len(agents_nonzero)/63336


#     agents_nonzero['tax_as_prop_of_pretax_profit'] = agents_nonzero.apply(lambda x: x['wq_pollutant_costs'] / (x['Yearly profit ($)']+x['wq_pollutant_costs']),axis=1)
#     agents_nonzero['tax_as_prop_of_pretax_profit'].median()
#     agents_nonzero['tax_as_prop_of_pretax_profit'] = agents_nonzero.apply(lambda x: x['wq_pollutant_costs'] / (max(0.001,x['Yearly profit ($)']+x['wq_pollutant_costs'])),axis=1)
#     agents_nonzero['tax_as_prop_of_pretax_profit'].median()
#     agents_nonzero['wq_pollutant_costs'].mean()
#     agents_nonzero['Yearly profit ($)'].mean()
#     (agents_nonzero['Yearly profit ($)']+agents_nonzero['wq_pollutant_costs']).mean()
    
#     prop = agents_nonzero.loc[(agents_nonzero['Yearly profit ($)']<0)&(agents_nonzero['wq_pollutant_costs']>(-1 * agents_nonzero['Yearly profit ($)']))]
#     len(prop)
#     len(prop)/len(agents_nonzero)
#     prop = agents_nonzero.loc[(agents_nonzero['Yearly profit ($)']<0)]
#     len(prop)
#     prp = grouped_to_agents.loc[(grouped_to_agents['Yearly profit ($)']<0)]
#     len(prp)
#     len(prp)/len(grouped_to_agents)


# initial = data[data["Timestep"]==0]
# initial['lu_class'] = initial["Land use category"].apply(lambda x: cat_to_class[x])
# initial['lucas_lu_class'] = initial["lucas_2016_land_use"].apply(lambda x: lucas_cat_to_class[x])
# sampled_initial = initial.sample(n=100000, random_state=42, replace=True)

# conf = confusion_matrix(sampled_initial['lu_class'], sampled_initial['lucas_lu_class'])
# conf_df = pd.DataFrame(conf)
# conf_df.columns = ["Other","Dairy","Sheep and beef","Horticulture","Arable"]
# conf_df.index = ["Other","Dairy","Sheep and beef","Horticulture","Arable"]
# conf_df.to_csv(f'{output_dir}plots/confusion_matrix.csv')

# data.loc[data["Land use"]=="carbon-forestry","Nitrogen loss (kg/ha/yr)"] = random.random() * 5
# relplot = sb.relplot(data=data[(data["Timestep"] == 25) & (data["Time in current land use (yr)"] > 1)],x="Profit ANPV ($/ha/yr)",y="Nitrogen loss (kg/ha/yr)",hue="Land use",kind="scatter")
# plt.savefig(f'{output_dir}plots/profit_vs_nloss.png')
# plt.close()

# first_and_last = data[(data.Timestep==1) | (data.Timestep==25)]
# for metric in ["Profit ANPV ($/ha/yr)","Nitrogen loss (kg/ha/yr)", "GHG Emissions (kg CO2e/ha/yr)","Nutrient Loss Mitigation","Intensity", "Utility ANPV ($/ha equivalent)"]:
#     plot = sb.displot(first_and_last,x=metric,hue="Land use",col="Timestep",bins=25)
#     plt.savefig(f'{output_dir}plots/{metric.split(" ")[0]}_hist.png')
#     plt.close()
#     no_bb_cherry = first_and_last[(first_and_last["Land use"] != "blueberry-covered") & (first_and_last["Land use"] != "cherry")]
#     plot = sb.displot(no_bb_cherry,x=metric,hue="Land use",col="Timestep",bins=25)
#     plt.savefig(f'{output_dir}plots/{metric.split(" ")[0]}_hist_no_bb_cherry.png')
#     plt.close()

# last = data[(data.Timestep==25)]
# first = data[(data.Timestep==1)]

# def create_pairplot(df,name,hue="Land use"):
#     pairplot = sb.pairplot(df,hue=hue)
#     plt.legend(title=None)
#     for i,ax in enumerate(pairplot.axes.flat):
#         # print(i)
#         yoffset = -0.2-0.025*(i%8)
#         xoffset = -0.1-0.1*(i%2)
#         ax.yaxis.set_label_coords(yoffset,0.5)
#         ax.xaxis.set_label_coords(0.5,xoffset)
#     plt.savefig(f'{output_dir}plots/{name}.png')
#     plt.close()

# learning_from_neighbours = sb.regplot(data=last,x="num_neighbours",y="Production Capability",lowess=True)
# plt.savefig(f'{output_dir}plots/learning_from_neighbours.png')
# plt.close()


# prof_df = last[["Land use","Production Capability","Intensity","Methane Emissions (kg CO2e/ha/yr)","Profit ANPV ($/ha/yr)"]]
# prof_df_init = first[["Land use","Production Capability","Intensity","Methane Emissions (kg CO2e/ha/yr)","Profit ANPV ($/ha/yr)"]]
# prof_df = prof_df[prof_df["Land use"] != "blueberry-covered"]
# create_pairplot(prof_df_init,"profit_pairplot_initial")
# create_pairplot(prof_df,"profit_pairplot_final")


# just_after_regs = data[data["Timestep"]==12]
# env_df = just_after_regs[["is_non_compliant","Intensity","nitrogen_weighting","Nutrient Loss Mitigation"]]
# create_pairplot(env_df,"env_pairplot_final",hue="is_non_compliant")

# ghg_df = last[["Land use","Intensity","Nitrous Oxide Emissions (kg CO2e/ha/yr)","Methane Emissions (kg CO2e/ha/yr)","Production Capability"]]
# create_pairplot(ghg_df,"ghg_pairplot_final")

# data["Relative utility"] = data.apply(lambda x: x["Utility ANPV ($/ha equivalent)"]*(1/max(data[data["Land use"]==x["Land use"]]["Utility ANPV ($/ha equivalent)"])),axis=1)
# data["likelihood_to_change"] = data['likelihood_to_change'].apply(lambda x: 0.5*round(2*x,1))
# grouped = data.groupby(["likelihood_to_change"]).mean().reset_index()
# relplot = sb.relplot(data=grouped,x="likelihood_to_change",y="Relative utility")
# plt.savefig(f'{output_dir}plots/likelihood_to_change.png')
# plt.close()


# cap_dist = pd.read_csv(f"outputs/{run_id}/cap_dist_1.csv")
# for lu in pd.unique(cap_dist["land_use_name"]):
#     filtered = cap_dist[(cap_dist["land_use_name"]==lu)]
#     plot = sb.displot(data=filtered,x="production_capability",hue="active",multiple='stack',bins=25)
#     plt.savefig(f'{output_dir}plots/cap_dist_{lu}.png')
#     plt.close()


# count = data.groupby(["Timestep","Land use"]).count().reset_index()
# # print(count)
# plot = sb.relplot(data=count,x="Timestep",y="uid",hue="Land use",kind="line")
# plt.savefig(f'{output_dir}plots/count.png')
# plt.close()


