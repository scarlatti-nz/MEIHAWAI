import pandas as pd
import scipy
from scipy.optimize import minimize

# fits coefficients for nutrient loss equations for each typology based on observed intensities in prior initialisation run. Parameters used in 'dairy typology' sheet of data/Model inputs.xlsx
# method deprecated for arable and horticultural land uses, opting for constant nutrient losses for each land use based on values in https://www.tandfonline.com/doi/full/10.1080/03036758.2022.2137532

dairy_typology_data = pd.read_excel("data/dairy-typology-table-for-n-p-ghg-losses-26-sept-2023.xlsx",sheet_name="LUT")
observed_intensities = pd.read_csv("observed_intensities.csv")

def objective_function(params,min,max,fifth,median,ninetyfifth):
    a,b,c = params
    loss_at_fifth = a + b * fifth + c * fifth**2
    loss_at_median = a + b * median + c * median**2
    loss_at_ninetyfifth = a + b * ninetyfifth + c * ninetyfifth**2
    obj = 0
    median_loss = (3 * min + max) / 4
    obj += (loss_at_fifth - min)**2
    obj += (loss_at_median - median_loss)**2
    obj += (loss_at_ninetyfifth - max)**2
    return obj

rows = []
for lu in pd.unique(observed_intensities['land_use_name']):
    for typology in range(1,25):
        if typology == 20:
            continue
        if lu == "dairy":
            N_scale = 1
            P_scale = 1
            N2O_scale = 1
        elif lu == "sheep-and-beef":
            N_scale = 0.25
            P_scale = 0.73
            N2O_scale = 0.2
        else:
            N_scale = 0.5
            P_scale = 0.5
            N2O_scale = 0.5
        
        min_N = dairy_typology_data[(dairy_typology_data['Typology number'] == typology)]['Min_N'].values[0] * N_scale
        max_N = dairy_typology_data[(dairy_typology_data['Typology number'] == typology)]['Max_N'].values[0] * N_scale
        min_P = dairy_typology_data[(dairy_typology_data['Typology number'] == typology)]['Min_P'].values[0] * P_scale
        max_P = dairy_typology_data[(dairy_typology_data['Typology number'] == typology)]['Max_P'].values[0] * P_scale
        min_N2O = dairy_typology_data[(dairy_typology_data['Typology number'] == typology)]['Min_N2O'].values[0] * N2O_scale
        max_N2O = dairy_typology_data[(dairy_typology_data['Typology number'] == typology)]['Max_N2O'].values[0] * N2O_scale
        fifth = observed_intensities[(observed_intensities['land_use_name'] == lu)]['fifth'].values[0]
        median = observed_intensities[(observed_intensities['land_use_name'] == lu)]['median'].values[0]
        ninetyfifth = observed_intensities[(observed_intensities['land_use_name'] == lu)]['ninetyfifth'].values[0]
        n_res = minimize(objective_function,[max_N,max_N*0.5,2],args=(min_N,max_N,fifth,median,ninetyfifth),bounds=[(0,max_N*2),(0,max_N),(0,max_N)])
        p_res = minimize(objective_function,[max_P,max_P*0.5,2],args=(min_P,max_P,fifth,median,ninetyfifth),bounds=[(0,max_P*2),(0,max_P),(0,max_P)])
        n2o_res = minimize(objective_function,[max_N2O,max_N2O*0.5,2],args=(min_N2O,max_N2O,fifth,median,ninetyfifth),bounds=[(0,max_N2O*2),(0,max_N2O),(0,max_N2O)])

        n_coefficients = n_res.x
        p_coefficients = p_res.x
        n2o_coefficients = n2o_res.x
        rows.append([lu,typology,fifth,ninetyfifth,n_coefficients[0],n_coefficients[1],n_coefficients[2],p_coefficients[0],p_coefficients[1],p_coefficients[2],n2o_coefficients[0],n2o_coefficients[1],n2o_coefficients[2]])
df = pd.DataFrame(rows,columns=['Land use','Typology','Intensity_lower_limit','Intensity_upper_limit','Constant_N','Linear_N','Quadratic_N','Constant_P','Linear_P','Quadratic_P','Constant_N2O','Linear_N2O','Quadratic_N2O'])
df.to_csv("fitted_nl_coefficients.csv")