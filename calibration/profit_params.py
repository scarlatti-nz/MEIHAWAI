import geopandas as gpd
import pandas as pd
import scipy as sp
from scipy import optimize
import numpy as np
import pyogrio
import time

# fits profit parameters for each land use based on observable revenue and profit found in max_revenue_profit_gdp.xlsx. Parameters used in 'land_use_parameters' sheet of data/Model inputs.xlsx

def calculate_observables(params,caps,observables,cap_x_transform):
    rev,varcost,h,fixedcost = params
    breakeven_cap,median_cap,maximal_cap = caps
    breakeven_cap = breakeven_cap + cap_x_transform
    median_cap = median_cap + cap_x_transform
    maximal_cap = maximal_cap + cap_x_transform
    minimizer = 0
    for cap,obs in zip([breakeven_cap,median_cap,maximal_cap],observables[:-1]):
        intensity = np.log(rev * cap * h / varcost) / (h * cap)
        penalty_for_intensity_over_limit = (max(0,intensity - 1) * obs)**3
        intensity = max(0,min(1,intensity))
        profit = rev * (1-np.e**(-h * cap * intensity)) - varcost * intensity - fixedcost
        minimizer += (profit-obs)**2
        minimizer += penalty_for_intensity_over_limit
    intensity = np.log(rev * median_cap * h / varcost) / (h * median_cap)
    intensity = max(0,min(1,intensity))
    rev_med = rev * (1-np.e**(-h * median_cap * intensity))
    minimizer += (rev_med - observables[-1])**2
    return minimizer
        

rev_data = pd.read_excel('max_revenue_profit_gdp.xlsx',sheet_name='max_revenue_profit_gdp')
crop_data = pd.read_csv("profit_params.csv")
dfrows=[]
for _,row in crop_data.iterrows():
    breakeven_cap = row["Breakeven capability"]
    median_cap = row["Median capability"]
    maximal_cap = row["Industry-leading capability"]
    ratio = row["Ratio of industry leading profit to median profit"]
    lu = row["Crop"]
    cap_x_transform = row['Capability x transform']
    # now fitting fixed cost rather than using known values
    # fixedcost = row["Fixed cost"]
    # annualisation of capital cost 
    breakeven_prof = float(row["Capital cost"]) / sum([1/1.06**i for i in range(1,51)])
    median_rev = rev_data.loc[rev_data['Crop']==lu,'Revenue'].values[0]
    median_prof = rev_data.loc[rev_data['Crop']==lu,'Profit'].values[0]
    maximal_prof =  ratio * median_prof
    observables = [breakeven_prof,median_prof,maximal_prof,median_rev]
    caps = [breakeven_cap,median_cap,maximal_cap]
    lbs = [median_rev,0.2*median_rev,1,0]
    ubs = [4*median_rev,3*median_rev,100,3*median_rev]
    bounds = optimize.Bounds(lbs,ubs)
    res = optimize.minimize(calculate_observables,x0=[1.2*median_rev,0.5*median_rev,5,0.5*median_rev],args=(caps,observables,cap_x_transform),bounds=bounds, tol=1e-40)
    rev,varcost,h,fixedcost = res.x
    err = res.fun / median_rev**2
    inflection_point = (varcost * np.e / rev) / h
    intensity_at_breakeven = np.log(rev * breakeven_cap * h / varcost) / (h * breakeven_cap)
    rev_adjustment = rev / median_rev
    dfrow = [lu,breakeven_cap,median_cap,maximal_cap,ratio,median_rev,median_prof,maximal_prof,rev,rev_adjustment,varcost,h,fixedcost,inflection_point,intensity_at_breakeven,err]
    dfrows.append(dfrow)
df = pd.DataFrame(dfrows,columns=["Crop","Breakeven capability","Median capability","Industry-leading capability","Ratio of industry leading profit to median profit","Median revenue","Median profit","Industry leader profit","Revenue coefficient","Revenue adjustment coefficient","Variable cost coefficient","Capability scalar","Fixed cost","Inflection point","Intensity at breakeven point","Error"])
while True:
    try:
        df.to_csv("profit_params_output.csv",index=False)
        break
    except:
        print("File open trying again")
        time.sleep(1)
        continue