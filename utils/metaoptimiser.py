import pandas as pd
from pandas.core.indexes import base
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
import os
import time

# baseline = "outputs/2024_02_16_18_00_00e0.0n0.0p0.0sed0.0n2o0.0ch40.0co270.0bfalsenolucfalsepathparabolic/plots/"
baseline = "../nesi_downloads/meihawai_outputs/n0p0ch40n2o0sed0escale0/plots/"

def objective_function(ch4,n2o,n,p,sed,baseline_profit,policy_profit,extension):
    obj = 0
    ghg_target = 0.65
    water_target = 0.1
    discount_rate = 0.06
    profit_scalar = 0.01 / baseline_profit[0]
    obj += max(0,(ch4-ghg_target)/(1-ghg_target))**2
    obj += max(0,(n2o-ghg_target)/(1-ghg_target))**2
    obj += max(0,(n-water_target)/(1-water_target))**2
    obj += max(0,(p-water_target)/(1-water_target))**2
    # obj += max(0,(sed-water_target)/(1-water_target))**2
    for i in range(len(baseline_profit)):
        obj += profit_scalar * ( baseline_profit[i] - policy_profit[i] + extension[i]) * 1/((1+discount_rate)**i)
    return obj
    
def get_obj_from_output(baseline,dir,profit_metric="producer"):
    
    outcomes = pd.read_csv(dir + "outcome_time_series.csv")
    
    ch4 = outcomes.loc[outcomes["Timestep"] == 25, "Methane Emissions"].iloc[0]
    n2o = outcomes.loc[outcomes["Timestep"] == 25, "Nitrous Oxide Emissions"].iloc[0]
    n = outcomes.loc[outcomes["Timestep"] == 25, "Out of compliance - N"].iloc[0]
    p = outcomes.loc[outcomes["Timestep"] == 25, "Out of compliance - P"].iloc[0]
    sed = outcomes.loc[outcomes["Timestep"] == 25, "Out of compliance - Sediment"].iloc[0]
    base_profit_csv = pd.read_csv(baseline + "outcomes_aggregated.csv")
    policy_profit_csv = pd.read_csv(dir + "outcomes_aggregated.csv")
    extension_csv = pd.read_csv(dir + "extension.csv")
    if profit_metric == "producer":
        baseline_profit = [base_profit_csv.loc[base_profit_csv["Timestep"] == i, "Yearly profit ($)"].sum() for i in range(25)]
        policy_profit = [policy_profit_csv.loc[policy_profit_csv["Timestep"] == i, "Yearly profit ($)"].sum() for i in range(25)]
    elif profit_metric == "welfare":
        base_profit_csv['Welfare ($)'] = base_profit_csv['Yearly profit ($)'] + base_profit_csv['wq_pollutant_costs'] + base_profit_csv['ghg_emission_costs'] - base_profit_csv['subsidy spent']
        policy_profit_csv['Welfare ($)'] = policy_profit_csv['Yearly profit ($)'] + policy_profit_csv['wq_pollutant_costs'] + policy_profit_csv['ghg_emission_costs'] - policy_profit_csv['subsidy spent']
        baseline_profit = [base_profit_csv.loc[base_profit_csv["Timestep"] == i, "Welfare ($)"].sum() for i in range(25)]
        policy_profit = [policy_profit_csv.loc[policy_profit_csv["Timestep"] == i, "Welfare ($)"].sum() for i in range(25)]
    extension_spend = [0.0] + [extension_csv.loc[(extension_csv["Timestep"] == i) & (extension_csv['funding']=='Additional') & (extension_csv['run_this_timestep'] == True), "cost"].sum() for i in range(1,25)]
    return objective_function(ch4,n2o,n,p,sed,baseline_profit,policy_profit,extension_spend)


def runner_func(prices):
    print(prices)
    n_price,p_price = prices
    ghg_price = 0.0
    sed_price = 0.0
    extension_scale = 0.0
    extension_scale = round(extension_scale,3)
    ch4price = round(ghg_price,3)
    n2oprice = round(ghg_price,3)
    n_price = round(n_price,3)
    p_price = round(p_price,3)
    sed_price = round(sed_price,3)
    timestamp_override = '2024_04_08_18_00_00'
    dir = "outputs/" + timestamp_override + "e" + str(extension_scale) + "n"+str(n_price)+"p"+str(p_price)+"sed"+str(sed_price)+"n2o"+str(n2oprice)+"ch4"+str(ch4price)+"co270.0bfalsenolucfalsepathparabolic/plots/"
    os.system("snakemake -c6 --config extensionScale=" + str(extension_scale) + " CH4Price=" + str(ch4price) + " N2OPrice=" + str(n2oprice) + " NPrice=" + str(n_price) + " PPrice=" + str(p_price) + " SedimentPrice=" + str(sed_price))
    obj = get_obj_from_output(baseline,dir)
    try:
        for i in range(26):
            os.remove(dir[:-6] + "agents_" + str(i) + ".arrow")
            os.remove(dir[:-6] + "catchments_" + str(i) + ".arrow")
            if i>0:
                os.remove(dir[:-6] + "extension_" + str(i) + ".arrow")
    except:
        pass
    return obj

# res = differential_evolution(runner_func, bounds=[(100,200),(500,2000)], popsize=4, maxiter=3, tol=0.01)

# print(res)

i = 0
rows = []
for dir in os.listdir("../nesi_downloads/meihawai_outputs/"):
# for dir in os.listdir('outputs/'):
    if ('2024' not in dir) and ('n2o' in dir):
        try:
            
            i += 1
            dir = "../nesi_downloads/meihawai_outputs/" + dir + "/plots/"
            obj = get_obj_from_output(baseline,dir)
            obj_welfare = get_obj_from_output(baseline,dir,profit_metric="welfare")
            base_profit_csv = pd.read_csv(baseline + "outcomes_aggregated.csv")
            policy_profit_csv = pd.read_csv(dir + "outcomes_aggregated.csv")
            baseline_profit = base_profit_csv["Yearly profit ($)"].sum()
            policy_profit = policy_profit_csv["Yearly profit ($)"].sum()
            profit_dip_abs = baseline_profit - policy_profit
            profit_dip_rel = profit_dip_abs / baseline_profit
            outcomes = pd.read_csv(dir + "outcome_time_series.csv")
            ch4 = outcomes.loc[outcomes["Timestep"] == 25, "Methane Emissions"].iloc[0]
            n2o = outcomes.loc[outcomes["Timestep"] == 25, "Nitrous Oxide Emissions"].iloc[0]
            n = outcomes.loc[outcomes["Timestep"] == 25, "Out of compliance - N"].iloc[0]
            p = outcomes.loc[outcomes["Timestep"] == 25, "Out of compliance - P"].iloc[0]
            sed = outcomes.loc[outcomes["Timestep"] == 25, "Out of compliance - Sediment"].iloc[0]
            nprice  = float(dir.split("/")[-3].split("n")[1].split("p")[0])
            pprice  = float(dir.split("/")[-3].split("p")[1].split("ch4")[0])
            ch4price = float(dir.split("/")[-3].split("ch4")[1].split("n2o")[0])
            n2oprice = float(dir.split("/")[-3].split("n2o")[1].split("sed")[0])
            sedprice = float(dir.split("/")[-3].split("sed")[1].split("escale")[0])
            escale = float(dir.split("/")[-3].split("escale")[1])
            # escale = float(dir.split("/")[-3].split("e")[1].split("n")[0])
            # escale = 0.0
            rows.append([nprice,pprice,ch4price,n2oprice,sedprice,n,p,n2o,ch4,sed,escale,profit_dip_abs,profit_dip_rel,obj,obj_welfare])
        except:
            pass
        
df = pd.DataFrame(rows, columns=["Nprice","Pprice","MethanePrice","N2O price","SedimentPrice","N target","P target","N2O target","CH4 target","Sediment target","Extension scale","25-year profit deficit (abs)","25-year profit deficit (rel)","Objective (producer profit)","Objective (welfare)"])
df.to_csv("outputs/nesiparamsearch3.csv",index=False)



