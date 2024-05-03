import matplotlib

matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['font.size'] = 14

import matplotlib.pyplot as plt



import seaborn as sb
import pandas as pd
import scipy.stats as stats
import random
import os
import argparse
import numpy as np
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import importlib
import math
import adjustText
import rasterio as rio
from rasterio.features import rasterize
import pyogrio
import scipy as sp
import pickle
import rasterio
import sys
# sb.set(font='Arial')
sb.set_palette("bright")


lu_to_cat = {"apple":"Arable/horticulture",
            "avocado":"Arable/horticulture",
            "blueberry-covered":"Arable/horticulture",
            "cherry":"Arable/horticulture",
            "kiwifruit-gold":"Arable/horticulture",
            "kiwifruit-green":"Arable/horticulture",
            "maize-grain":"Arable/horticulture",
            "onions":"Arable/horticulture",
            "peas-vining":"Arable/horticulture",
            "potatoes":"Arable/horticulture",
            "wheat":"Arable/horticulture",
            "wine-grape-pinot-noir":"Arable/horticulture",
            "wine-grape-sauvignon-blanc":"Arable/horticulture",
            "dairy":"Dairy",
            "dairy-irrigated":"Dairy",
            "sheep-and-beef":"Sheep and beef",
            "dairy-bolus":"Dairy",
            "sheep-and-beef-bolus":"Sheep and beef",
            "production-forestry":"Forestry",
            "carbon-forestry":"Forestry",
            "native-forestry":"Forestry",
            "fallow":"Fallow/bare land"}


output_dir = sys.argv[1]
baseline_dir = output_dir + "/status_quo/plots/"
fully_optimised_dir = output_dir + "/policy_optimised/plots/"
zerocarbon = output_dir + "/no_carbon_price/plots/"
extensiononlylow = output_dir + "/status_quo_plus_extension/plots/"
figures_dir = sys.argv[2]
    
# figures_dir = "plots/pricing/"

try:
    os.mkdir(figures_dir)
except:
    pass
run_scale = 1.0

# numbers_for_quarto = {

#     "Nitrogen price" : fully_optimised_dir.split("/")[-3].split("n")[1].split("p")[0],
#     "Phosphorous price" : fully_optimised_dir.split("/")[-3].split("p")[1].split("sed")[0],
#     # "Sediment price" : fully_optimised_dir.split("sed")[1].split("n2o")[0],
#     "Nitrous oxide price" : fully_optimised_dir.split("/")[-3].split("n2o")[1].split("ch4")[0],
#     "Methane price" : fully_optimised_dir.split("/")[-3].split("ch4")[1].split("co2")[0],
# }


def get_data(dir):
    outcomes = pd.read_csv(dir+'outcome_time_series.csv').drop(columns=["Unnamed: 0"])
    outcomes["Timestep"] += 2025
    outcomes.rename(columns={"Methane Emissions":"Methane","Nitrous Oxide Emissions":"Nitrous oxide","Out of compliance - N":"Nitrogen","Out of compliance - P":"Phosphorous","Out of compliance - Sediment":"Sediment"},inplace=True)
    # removing sediment for aares figs
    outcomes.drop("Sediment",axis=1,inplace=True)
    outcomes[[c for c in outcomes.columns if "Timestep" != c]] *= 100

    aggregates = pd.read_csv(dir+'outcomes_aggregated.csv').drop(columns=["Unnamed: 0"])
    aggregates = aggregates[aggregates["Land use"] != "fallow"]
    aggregates["Timestep"] += 2025
    if "switching_costs_incurred" in aggregates.columns:
        for lu in aggregates["Land use"].unique():
            for timestep in range(2025,2051):
                total_sc_depreciated = aggregates.loc[(aggregates["Timestep"] == timestep) & (aggregates["Land use"] == lu), "switching_costs_incurred"].sum() / 20
                aggregates.loc[(aggregates["Land use"] == lu) & (aggregates["Timestep"] >= timestep) & (aggregates["Timestep"] <= timestep+20), "Yearly profit ($)"] -= total_sc_depreciated
    else:
        print("WARNING NOT CALCULATING SWITCHING COSTS DEPRECIATION")
    aggregates["Net Welfare"] = aggregates["Yearly profit ($)"] + aggregates['ghg_emission_costs'] + aggregates["wq_pollutant_costs"] - aggregates['subsidy spent']

    
    aggregates["Land use category"] = aggregates["Land use"].map(lu_to_cat)
    aggregates_by_cat = aggregates.groupby(["Land use category","Timestep"]).sum().reset_index()
    forestry_zero_row = {col: 0 for  col in aggregates_by_cat.columns}
    forestry_zero_row["Land use category"] = "Forestry"
    forestry_zero_row["Timestep"] = 2025
    aggregates_by_cat = aggregates_by_cat.append(forestry_zero_row,ignore_index=True)
    # adjust for existing forestry area, scaled to model run scale
    aggregates_by_cat.loc[aggregates_by_cat["Land use category"] == "Forestry","Area (ha)"] += 1.6e6
    aggregates_by_cat.sort_values(by=["Land use category","Timestep"],inplace=True)

    return outcomes, aggregates, aggregates_by_cat


def assign_endlabels(plot,df,xlabel,ylabel,hue,ylim,lpad=1):
    handles, labels = plot.get_legend_handles_labels()
    texts = []
    endpoints = sorted(df.loc[df[xlabel] == df[xlabel].max(),ylabel].values)
    adjusted = endpoints.copy()
    minygaps = [0.12 * ylim if '\n' in label else 0.06 * ylim for label in labels]
    # if any('\n' in label for label in labels):
    #     minygap = ylim * 0.12
    # else:
    #     minygap = ylim * 0.06
    while any([adjusted[i+1] - adjusted[i] < minygaps[i+1] for i in range(len(adjusted)-1)]):
        for i in range(len(adjusted)-1):
            minygap = minygaps[i+1]
            diff = adjusted[i+1] - adjusted[i]
            if diff < minygap:
                adjusted[i+1] = adjusted[i+1] + minygap/2
                adjusted[i] = adjusted[i] - minygap/2
        
        adjusted[0] -= 1
    yoffsetdict = dict(zip(endpoints,adjusted))


    for handle, label in zip(handles, labels):
        last_valid = df[df[hue] == label].dropna().iloc[-1]
        y = last_valid[ylabel]
        x = last_valid[xlabel]
        t = plot.text(x, yoffsetdict[y], ' '*lpad + label, color=handle.get_color(), 
                verticalalignment='center', bbox=dict(facecolor='white', alpha=0.0, edgecolor='none', boxstyle='round,pad=0.1'))
        texts.append(t)
    plot.legend().remove()



def fig1(base,comp):
    baseline_outcomes, baseline_aggregates, baseline_aggregates_by_cat = get_data(base)
    comparison_outcomes, comparison_aggregates, comparison_aggregates_by_cat = get_data(comp)
    fig, axs = plt.subplots(1,2, figsize=(10, 5), gridspec_kw={'wspace':0.8, 'hspace':0.3, 'right':0.8, 'left':0.1})
    baseline_aggregates_by_cat["Area (000 ha)"] = baseline_aggregates_by_cat["Area (ha)"] / 1e3
    comparison_aggregates_by_cat["Area (000 ha)"] = comparison_aggregates_by_cat["Area (ha)"] / 1e3
    
    ylim = math.ceil(max(baseline_aggregates_by_cat["Area (000 ha)"].max(),comparison_aggregates_by_cat["Area (000 ha)"].max())/1e2) * 1e2
    for i in range(2):
        ax = axs[i]
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.text(0, 1.15, chr(97 + i), transform=ax.transAxes, 
            fontsize=18, fontweight='bold', va='top', ha='right')
        ax.set_ylim(0,ylim)
        ax.set_xlim(2025,2050)

    ax = axs[0]
    plot = sb.lineplot(data=baseline_aggregates_by_cat,x="Timestep",y="Area (000 ha)",hue="Land use category",ax=ax)
    assign_endlabels(plot,baseline_aggregates_by_cat,"Timestep","Area (000 ha)","Land use category",ylim,lpad=1)
    ax.set_xlabel("")
    ax.title.set_text("Land use (Status quo)")
    ax.title.set_ha("left")
    ax.title.set_position([0.0,1.0])

    ax = axs[1]
    plot = sb.lineplot(data=comparison_aggregates_by_cat,x="Timestep",y="Area (000 ha)",hue="Land use category",ax=ax)
    assign_endlabels(plot,comparison_aggregates_by_cat,"Timestep","Area (000 ha)","Land use category",ylim,lpad=1)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_yticklabels([])
    ax.title.set_text("Land use (Extension & incentives)")
    ax.title.set_ha("left")
    ax.title.set_position([0.0,1.0])
    plt.savefig(f'{figures_dir}/results_fig_1.png',dpi=300)
    plt.close()

    # isolate fig 1a
    plot = sb.lineplot(data=baseline_aggregates_by_cat,x="Timestep",y="Area (000 ha)",hue="Land use category")
    assign_endlabels(plot,baseline_aggregates_by_cat,"Timestep","Area (000 ha)","Land use category",ylim,lpad=1)
    plt.subplots_adjust(right=0.7)
    plot.set_xlabel("")
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plot.set_ylim(0,ylim)
    plot.set_xlim(2025,2050)
    plt.savefig(f'{figures_dir}/results_fig_1a.png',dpi=300)
    plt.close()


def raster_comparison(start,end):
    if start == -2:
        return -1
    # how should we represent land not modelled?
    elif start == 0:
        return 0
    elif start == end:
        return 0
    elif end == -1:
        return 1
    else:
        return 2

def get_diff_of_diffs(baseline,comparison):
    if baseline == -1 or comparison == -1:
        return -1
    elif baseline == comparison:
        return 0
    elif comparison == 1:
        return 1
    else:
        return 2
    

def fig2v2(base):
    with rio.open(f'data/lucas_raster.tif') as r:
            lucas_raster = r.read(1)
    with rio.open(f'{base}initial_land_use.tif') as r:
            baseline_init = r.read(1)
    fig, axs = plt.subplots(1,2, figsize=(10, 8))   

    ax = axs[0]
    plot = sb.heatmap(lucas_raster,cmap=["white","lightgrey","lightgrey","orange","green","blue","red"],cbar=False,ax=ax)
    ax.title.set_text("Real world land use (LUCAS 2016)")
    ax.title.set_ha("left")
    ax.title.set_position([-0.05,1.0])

    ax = axs[1]
    plot = sb.heatmap(baseline_init,cmap=["white","lightgrey","lightgrey","orange","green","blue","red"],cbar=False,ax=ax)
    ax.title.set_text("Modelled land use (MEIHAWAI initialisation)")
    ax.title.set_ha("left")
    ax.title.set_position([-0.05,1.0])

    for i in range(2):
        ax = axs[i]
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.text(0, 1.1, chr(97 + i), transform=ax.transAxes, 
            fontsize=18, fontweight='bold', va='top', ha='right')
        ax.set_xticklabels([])
        ax.set_xticks([])
        ax.set_yticklabels([])
        ax.set_yticks([])
        # ax.set_ylim(900,0)
    plt.savefig(f'{figures_dir}/lucas_v_meihawai.png',dpi=300)
    plt.close()

def fig3(base,comp,suffix=""):
    baseline_outcomes, baseline_aggregates, baseline_aggregates_by_cat = get_data(base)
    comparison_outcomes, comparison_aggregates, comparison_aggregates_by_cat = get_data(comp)
    sqname = "Status quo"
    if suffix == "":
        scenario_name = "Extension & incentives"
    elif suffix == "_zerocarbon":
        scenario_name = "No carbon price"
    elif suffix == "_extensiononlylow":
        scenario_name = "Extension only"
    elif suffix == "_extensiononlyhigh":
        scenario_name = "Extension only"
    elif suffix == "_pricingandextension":
        scenario_name = "Extension & incentives"
        sqname = "Incentives only"

    fig, axs = plt.subplots(2, 2, figsize=(10, 10), gridspec_kw={'wspace':1.2, 'hspace':0.3, 'right':0.75, 'left':0.12})
    
    for i in range(4):
        ax = axs[i//2, i%2]
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.text(0, 1.2, chr(97 + i), transform=ax.transAxes, 
            fontsize=18, fontweight='bold', va='top', ha='right')
        ax.set_xlim(2025,2050)

    ax = axs[0,0]
    ax.set_ylim(0,120)
    melted = pd.melt(baseline_outcomes,id_vars=["Timestep"])
    plot = sb.lineplot(data=melted,x="Timestep",y="value",hue="variable",ax=ax)
    ax.axhline(y=65,linestyle=':',color='black')
    ax.text(2035, 60, r'CH$_4$, N$_2$O target', horizontalalignment='center', verticalalignment='center')
    ax.axhline(y=10,linestyle=':',color='black')
    ax.text(2037, 15, 'N, P, Sediment target', horizontalalignment='center', verticalalignment='center')
    assign_endlabels(plot,melted,"Timestep","value","variable",120)
    ax.set_xlabel("")
    ax.set_ylabel("Pollution level index")
    ax.set_xticklabels([])
    ax.title.set_text(f"Pollution ({sqname})")
    ax.title.set_ha("left")
    ax.title.set_position([0.0,1.0])

    ax = axs[0,1]
    ax.set_ylim(0,120)
    melted = pd.melt(comparison_outcomes,id_vars=["Timestep"])
    plot = sb.lineplot(data=melted,x="Timestep",y="value",hue="variable",ax=ax)
    ax.axhline(y=65,linestyle=':',color='black')
    # ax.text(2035, 55, 'CH4, N2O target', horizontalalignment='center', verticalalignment='center')
    ax.axhline(y=10,linestyle=':',color='black')
    # ax.text(2035, 20, 'N, P, Sediment target', horizontalalignment='center', verticalalignment='center')
    assign_endlabels(plot,melted,"Timestep","value","variable",120)
    ax.set_xlabel("")
    ax.set_ylabel("Pollution level index")
    ax.set_xticklabels([])
    ax.set_yticklabels([])      
    ax.title.set_text(f"Pollution ({scenario_name})")
    ax.title.set_ha("left")
    ax.title.set_position([0.0,1.0])

    ax = axs[1,0]
    land_use_comparison = baseline_aggregates_by_cat.merge(comparison_aggregates_by_cat,on=["Land use category","Timestep"],suffixes=["_baseline","_comparison"])
    # land_use_comparison[f"Land use index\n({sqname} scenario = 100)"] =  100 * land_use_comparison["Area (ha)_comparison"] / land_use_comparison["Area (ha)_baseline"]
    land_use_comparison[f"Change in land use (000 ha)\nrelative to {sqname} scenario"] =  1e-3 * (land_use_comparison["Area (ha)_comparison"] - land_use_comparison["Area (ha)_baseline"])
    plot = sb.lineplot(data=land_use_comparison,x="Timestep",y=f"Change in land use (000 ha)\nrelative to {sqname} scenario",hue="Land use category",ax=ax)
    ax.axhline(y=0,linestyle=':',color='black')
    # if suffix == "_zerocarbon":
    #     ax.set_ylim(-3000,3000)
    #     ygap = 6000
    # elif "extension" in suffix:
    #     ax.set_ylim(-50,50)
    #     ygap = 80
    # else:
    #     ax.set_ylim(-700,700)
    #     ygap = 800
    yheight = max(abs(land_use_comparison[f"Change in land use (000 ha)\nrelative to {sqname} scenario"]))
    print("yheight: ",yheight)
    ylim = round(yheight * 1.1, -1 * math.floor(math.log10(yheight)))
    ax.set_ylim(-1*ylim,ylim)
    ygap = 2*ylim
    assign_endlabels(plot,land_use_comparison,"Timestep",f"Change in land use (000 ha)\nrelative to {sqname} scenario","Land use category",ygap)
    ax.set_xlabel("")
    ax.title.set_text(f"Relative land use\n{scenario_name} vs {sqname}")
    ax.title.set_ha("left")
    ax.title.set_position([0.0,1.0])

    ax = axs[1,1]
    aggregate_profit = land_use_comparison.groupby("Timestep").sum().reset_index()
    aggregate_profit[f"{sqname} (profit)"] = 100 * aggregate_profit["Yearly profit ($)_baseline"] / aggregate_profit.loc[aggregate_profit["Timestep"]==2025,"Yearly profit ($)_baseline"].iloc[0]
    aggregate_profit[f"{scenario_name} (profit)"] = 100 * aggregate_profit["Yearly profit ($)_comparison"] / aggregate_profit.loc[aggregate_profit["Timestep"]==2025,"Yearly profit ($)_baseline"].iloc[0]
    print("profit index: ",aggregate_profit.loc[aggregate_profit["Timestep"]==2026,f"{scenario_name} (profit)"])
    
    
    # numbers_for_quarto["Extension spend"] = total_extension_spend
    # numbers_for_quarto["Public extension spend"] = public_extension_spend
    extension_csv = pd.read_csv(f'{comp}extension.csv')
    extension_spend = [0.0] + [extension_csv.loc[(extension_csv["Timestep"] == i) & (extension_csv['funding']=='Additional') & (extension_csv['run_this_timestep'] == True), "cost"].sum() for i in range(1,26)]
    print("extension spend: ", extension_spend)
    aggregate_profit["Net Welfare_comparison"] -= extension_spend
    
    aggregate_profit[f"{scenario_name} \n (profit + net tax)\n"] = 100 * aggregate_profit["Net Welfare_comparison"] / aggregate_profit.loc[aggregate_profit["Timestep"]==2025,"Yearly profit ($)_baseline"].iloc[0]
    if "pricing" in suffix:
        aggregate_profit[f"{sqname} \n (profit + net tax)\n"] = 100 * aggregate_profit["Net Welfare_baseline"] / aggregate_profit.loc[aggregate_profit["Timestep"]==2025,"Yearly profit ($)_baseline"].iloc[0]
        aggregate_profit.rename(columns={f"{scenario_name} (profit)":f"{scenario_name} \n (profit)"},inplace=True)
        melted = pd.melt(aggregate_profit,id_vars=["Timestep"],value_vars=[f"{sqname} (profit)",f"{scenario_name} \n (profit)", f"{sqname} \n (profit + net tax)\n" , f"{scenario_name} \n (profit + net tax)\n"])
    elif suffix == "":
        aggregate_profit.rename(columns={f"{scenario_name} (profit)":f"{scenario_name} \n (profit)"},inplace=True)
        melted = pd.melt(aggregate_profit,id_vars=["Timestep"],value_vars=[f"{sqname} (profit)",f"{scenario_name} \n (profit)", f"{scenario_name} \n (profit + net tax)\n"])
    else:
        melted = pd.melt(aggregate_profit,id_vars=["Timestep"],value_vars=[f"{sqname} (profit)",f"{scenario_name} (profit)", f"{scenario_name} \n (profit + net tax)\n"])
    print("Producer delta: ", aggregate_profit["Yearly profit ($)_baseline"].sum() - aggregate_profit["Yearly profit ($)_comparison"].sum())
    print("Public delta: ", aggregate_profit["Net Welfare_baseline"].sum() - aggregate_profit["Net Welfare_comparison"].sum())
    
    if suffix == "_zerocarbon":
        melted = melted[melted["variable"] != f"{scenario_name} \n (profit + net tax)\n"]
    
    # print(melted)
    plot = sb.lineplot(data=melted,x="Timestep",y="value",hue="variable",ax=ax,palette=sb.color_palette()[7:] + [sb.color_palette()[6]])
    ax.set_ylim(0,400)
    assign_endlabels(plot,melted,"Timestep","value","variable",400)
    ax.set_xlabel("")
    ax.set_ylabel("Aggregate welfare index")
    ax.title.set_text("Aggregate welfare indicators")
    ax.title.set_ha("left")
    ax.title.set_position([0.0,1.0])
    
    plt.savefig(f'{figures_dir}/results_fig_3{suffix}.png',dpi=300)
    plt.close()


    melted = pd.melt(baseline_outcomes,id_vars=["Timestep"])
    plot = sb.lineplot(data=melted,x="Timestep",y="value",hue="variable")
    plt.subplots_adjust(right=0.7)
    plot.axhline(y=65,linestyle=':',color='black')
    plot.text(2035, 60, r'CH$_4$, N$_2$O target', fontname='Arial', horizontalalignment='center', verticalalignment='center')
    plot.axhline(y=10,linestyle=':',color='black')
    plot.text(2037, 15, 'N, P, Sediment target', horizontalalignment='center', verticalalignment='center')
    assign_endlabels(plot,melted,"Timestep","value","variable",120)
    plot.set_xlabel("")
    plot.set_ylabel("Pollution level index")
    plot.set_xlim(2025,2050)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plt.savefig(f'{figures_dir}/results_fig_3a.png',dpi=300)
    plt.close()

    # if not suffix:
    #     plot = sb.lineplot(data=land_use_comparison,x="Timestep",y="Land use index\n(status quo scenario = 100)",hue="Land use category")
    #     plot.axhline(y=100,linestyle=':',color='black')
    #     plot.set_ylim(90,110)
    #     assign_endlabels(plot,land_use_comparison,"Timestep","Land use index\n(status quo scenario = 100)","Land use category",20)
    #     plot.set_xlabel("")
    #     # plot.set_ylabel("Relative land use index")
    #     plt.subplots_adjust(right=0.7,left=0.2)
    #     plot.set_xlim(2025,2050)
    #     plot.spines['right'].set_visible(False)
    #     plot.spines['top'].set_visible(False)
    #     plt.savefig(f'{figures_dir}/results_fig_3c.png',dpi=300)
    #     plt.close()

    #     melted = pd.melt(aggregate_profit,id_vars=["Timestep"],value_vars=["Status quo (profit)","Policy (profit)", "Policy \n (profit + net tax)\n"])
    #     plot = sb.lineplot(data=melted,x="Timestep",y="value",hue="variable",palette=sb.color_palette()[7:])
    #     plot.set_ylim(0,400)
    #     assign_endlabels(plot,melted,"Timestep","value","variable",400)
    #     plot.set_xlabel("")
    #     plot.set_ylabel("Aggregate welfare index")
    #     plot.spines['right'].set_visible(False)
    #     plot.spines['top'].set_visible(False)
    #     plt.subplots_adjust(right=0.7)
    #     plt.savefig(f'{figures_dir}/results_fig_3d.png',dpi=300)
    #     plt.close()




def fig3aares(zerocarbon,baseline,pricing,pricingandextension):
    _,_,zerocarbon_aggregates_by_cat = get_data(zerocarbon)
    _,_,baseline_aggregates_by_cat = get_data(baseline)
    _,_,pricing_aggregates_by_cat = get_data(pricing)
    _,_,pricingandextension_aggregates_by_cat = get_data(pricingandextension)
    
    land_use_comparison = baseline_aggregates_by_cat.merge(pricing_aggregates_by_cat,on=["Land use category","Timestep"],suffixes=["_baseline","_pricing"])
    zerocarbon_aggregates_by_cat.rename(lambda x: x+"_zerocarbon",axis=1,inplace=True)
    zerocarbon_aggregates_by_cat.rename(columns={"Land use category_zerocarbon":"Land use category","Timestep_zerocarbon":"Timestep"},inplace=True)
    pricingandextension_aggregates_by_cat.rename(lambda x: x+"_pricingandextension",axis=1,inplace=True)
    pricingandextension_aggregates_by_cat.rename(columns={"Land use category_pricingandextension":"Land use category","Timestep_pricingandextension":"Timestep"},inplace=True)
    
    land_use_comparison = land_use_comparison.merge(zerocarbon_aggregates_by_cat,on=["Land use category","Timestep"])
    land_use_comparison = land_use_comparison.merge(pricingandextension_aggregates_by_cat,on=["Land use category","Timestep"])


    palette = [sb.color_palette()[2]] + sb.color_palette()[7:]

    aggregate_profit = land_use_comparison.groupby("Timestep").sum().reset_index()
    aggregate_profit["No carbon price"] = 100 * aggregate_profit["Yearly profit ($)_zerocarbon"] / aggregate_profit.loc[aggregate_profit["Timestep"]==2025,"Yearly profit ($)_baseline"].iloc[0]
    aggregate_profit["Status quo"] = 100 * aggregate_profit["Yearly profit ($)_baseline"] / aggregate_profit.loc[aggregate_profit["Timestep"]==2025,"Yearly profit ($)_baseline"].iloc[0]
    aggregate_profit["Pricing only"] = 100 * aggregate_profit["Yearly profit ($)_pricing"] / aggregate_profit.loc[aggregate_profit["Timestep"]==2025,"Yearly profit ($)_baseline"].iloc[0]
    aggregate_profit["Pricing plus extension"] = 100 * aggregate_profit["Yearly profit ($)_pricingandextension"] / aggregate_profit.loc[aggregate_profit["Timestep"]==2025,"Yearly profit ($)_baseline"].iloc[0]
    melted = pd.melt(aggregate_profit,id_vars=["Timestep"],value_vars=["No carbon price", "Status quo","Pricing only", "Pricing plus extension"])
    plot = sb.lineplot(data=melted,x="Timestep",y="value",hue="variable",palette=palette)
    assign_endlabels(plot,melted,"Timestep","value","variable",400)
    plot.set_xlabel("")
    plot.set_ylabel("Aggregate producer profit")
    plot.set_xlim(2025,2050)
    plot.set_ylim(0,400)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plt.subplots_adjust(right=0.65)
    plt.savefig(f'{figures_dir}/results_fig_3b.png',dpi=300)
    plt.close()

    base_total_extension_spend  = 2e8 * float(baseline.split('/')[-3].split("e")[1].split("n")[0]) * run_scale
    base_public_extension_spend = get_public_spend(base_total_extension_spend,1.3e8 * run_scale)
    pricingandextension_total_extension_spend  = 2e8 * float(pricingandextension.split('/')[-3].split("e")[1].split("n")[0]) * run_scale
    pricingandextension_public_extension_spend = get_public_spend(pricingandextension_total_extension_spend,1.3e8 * run_scale)
    
    aggregate_profit["Net Welfare_zerocarbon"] -= base_public_extension_spend
    aggregate_profit["Net Welfare_baseline"] -= base_public_extension_spend
    aggregate_profit["Net Welfare_pricing"] -= base_public_extension_spend
    aggregate_profit["Net Welfare_pricingandextension"] -= pricingandextension_public_extension_spend

    aggregate_welfare = aggregate_profit.copy()

    aggregate_welfare["No carbon price"] = 100 * aggregate_welfare["Net Welfare_zerocarbon"] / aggregate_welfare.loc[aggregate_welfare["Timestep"]==2025,"Yearly profit ($)_baseline"].iloc[0]
    aggregate_welfare["Status quo"] = 100 * aggregate_welfare["Net Welfare_baseline"] / aggregate_welfare.loc[aggregate_welfare["Timestep"]==2025,"Yearly profit ($)_baseline"].iloc[0]
    aggregate_welfare["Pricing only"] = 100 * aggregate_welfare["Net Welfare_pricing"] / aggregate_welfare.loc[aggregate_welfare["Timestep"]==2025,"Yearly profit ($)_baseline"].iloc[0]
    aggregate_welfare["Pricing plus extension"] = 100 * aggregate_welfare["Net Welfare_pricingandextension"] / aggregate_welfare.loc[aggregate_welfare["Timestep"]==2025,"Yearly profit ($)_baseline"].iloc[0]
    melted = pd.melt(aggregate_welfare,id_vars=["Timestep"],value_vars=["No carbon price", "Status quo","Pricing only", "Pricing plus extension"])
    
    plot = sb.lineplot(data=melted,x="Timestep",y="value",hue="variable",palette=palette)
    assign_endlabels(plot,melted,"Timestep","value","variable",400)
    plot.set_xlabel("")
    plot.set_ylabel("Aggregate profit + net tax")
    plot.set_xlim(2025,2050)
    plot.set_ylim(0,400)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plt.subplots_adjust(right=0.65)
    # print(melted)
    plt.savefig(f'{figures_dir}/results_fig_3c.png',dpi=300)
    plt.close()

def fig4aares(pande):
    _,_,aggs = get_data(pande)
    initial_and_final = aggs.loc[(aggs["Timestep"]==2025) | (aggs["Timestep"]==2050)]
    initial_and_final["Land use category"] = initial_and_final["Land use category"].apply(lambda x: "New forestry" if x == "Forestry" else x)
    initial_and_final['Yearly profit ($)'] = initial_and_final['Yearly profit ($)'] / 1e9
    # Plot
    plot = sb.barplot(data=initial_and_final,y="Timestep",x="Yearly profit ($)",hue="Land use category",orient="h")
    plot.set_xlabel("Aggregate producer profit ($ Billion)")
    plot.set_ylabel("")
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)

    # plot = sb.objects.Plot(initial_and_final["Timestep"],initial_and_final["Yearly profit ($)"],initial_and_final["Land use category"]).add(sb.objects.Bar())

    plt.savefig(f'{figures_dir}/results_fig_4a.png',dpi=300)
    plt.close()


def create_deltas_df(dir):
    deltas = pd.read_csv(dir+'changed_lu.csv').drop(columns=["Unnamed: 0"])
    delta_N_changed_lu = deltas.loc[(deltas["N_catchment"]==True)&(deltas["stayed_in_same_lu"]==False)&(deltas["Timestep"]==25),"Total Nitrogen loss (kg/yr)"].sum() - deltas.loc[(deltas["N_catchment"]==True)&(deltas["stayed_in_same_lu"]==False)&(deltas["Timestep"]==0),"Total Nitrogen loss (kg/yr)"].sum()
    delta_N_same_lu = deltas.loc[(deltas["N_catchment"]==True)&(deltas["stayed_in_same_lu"]==True)&(deltas["Timestep"]==25),"Total Nitrogen loss (kg/yr)"].sum() - deltas.loc[(deltas["N_catchment"]==True)&(deltas["stayed_in_same_lu"]==True)&(deltas["Timestep"]==0),"Total Nitrogen loss (kg/yr)"].sum()
    delta_P_changed_lu = deltas.loc[(deltas["P_catchment"]==True)&(deltas["stayed_in_same_lu"]==False)&(deltas["Timestep"]==25),"Total Phosphorous loss (kg/yr)"].sum() - deltas.loc[(deltas["P_catchment"]==True)&(deltas["stayed_in_same_lu"]==False)&(deltas["Timestep"]==0),"Total Phosphorous loss (kg/yr)"].sum()
    delta_P_same_lu = deltas.loc[(deltas["P_catchment"]==True)&(deltas["stayed_in_same_lu"]==True)&(deltas["Timestep"]==25),"Total Phosphorous loss (kg/yr)"].sum() - deltas.loc[(deltas["P_catchment"]==True)&(deltas["stayed_in_same_lu"]==True)&(deltas["Timestep"]==0),"Total Phosphorous loss (kg/yr)"].sum()
    delta_methane_changed_lu = deltas.loc[(deltas["stayed_in_same_lu"]==False)&(deltas["Timestep"]==25),"Total Methane Emissions (kg CO2e/yr)"].sum() - deltas.loc[(deltas["stayed_in_same_lu"]==False)&(deltas["Timestep"]==0),"Total Methane Emissions (kg CO2e/yr)"].sum()
    delta_methane_same_lu = deltas.loc[(deltas["stayed_in_same_lu"]==True)&(deltas["Timestep"]==25),"Total Methane Emissions (kg CO2e/yr)"].sum() - deltas.loc[(deltas["stayed_in_same_lu"]==True)&(deltas["Timestep"]==0),"Total Methane Emissions (kg CO2e/yr)"].sum()
    delta_nitrous_changed_lu = deltas.loc[(deltas["stayed_in_same_lu"]==False)&(deltas["Timestep"]==25),"Total Nitrous Oxide Emissions (kg CO2e/yr)"].sum() - deltas.loc[(deltas["stayed_in_same_lu"]==False)&(deltas["Timestep"]==0),"Total Nitrous Oxide Emissions (kg CO2e/yr)"].sum()
    delta_nitrous_same_lu = deltas.loc[(deltas["stayed_in_same_lu"]==True)&(deltas["Timestep"]==25),"Total Nitrous Oxide Emissions (kg CO2e/yr)"].sum() - deltas.loc[(deltas["stayed_in_same_lu"]==True)&(deltas["Timestep"]==0),"Total Nitrous Oxide Emissions (kg CO2e/yr)"].sum()
    area_changed_lu = deltas.loc[(deltas["stayed_in_same_lu"]==False)&(deltas["Timestep"]==25),"Area (ha)"].sum()
    area_same_lu  = deltas.loc[(deltas["stayed_in_same_lu"]==True)&(deltas["Timestep"]==25),"Area (ha)"].sum()
    prop_changed = area_changed_lu / (area_changed_lu + area_same_lu)

    df_for_plotting = pd.DataFrame({"Pollutant":["N","P",r"CH$_4$",r"N$_2$O"],"Land use change":[delta_N_changed_lu,delta_P_changed_lu,delta_methane_changed_lu,delta_nitrous_changed_lu],"Management change":[delta_N_same_lu,delta_P_same_lu,delta_methane_same_lu,delta_nitrous_same_lu]})
    df_for_plotting[["Land use change","Management change"]] = df_for_plotting[["Land use change","Management change"]].div(df_for_plotting.sum(axis=1), axis=0)
    
    return df_for_plotting, prop_changed

def fig4(base,comp):
    fig, axs = plt.subplots(1,2,figsize=(8,5), gridspec_kw={'wspace':0.3, 'hspace':0.5, 'right':0.9, 'left':0.15, 'bottom':0.2})
    baseline, base_prop_changed_lu = create_deltas_df(base)
    policy, comp_prop_changed_lu = create_deltas_df(comp)
    print(base_prop_changed_lu,comp_prop_changed_lu)
    # numbers_for_quarto['Proportion of area switching land use (status quo)'] = base_prop_changed_lu
    # numbers_for_quarto['Proportion of area switching land use (policy)'] = comp_prop_changed_lu

    ax = axs[0]
    baseline.plot(kind='bar', stacked=False, ax=ax)
    ax.set_ylabel("Proportion of load reduction")
    ax.legend().remove()
    ax.set_xticklabels(policy["Pollutant"],rotation=0)
    ax.title.set_text("Status quo")
    ax.title.set_ha("left")
    ax.title.set_position([0.0,1.0])

    ax = axs[1]
    policy.plot(kind='bar', stacked=False, ax=ax)
    ax.set_xticklabels(policy["Pollutant"],rotation=0)
    ax.set_yticklabels([])
    ax.legend(bbox_to_anchor=(0.8,-0.2), loc='right', borderaxespad=0.1,ncol=2)
    # ax.set_ylabel("Proportion of load reduction")
    
    ax.title.set_text("Extension & incentives")
    ax.title.set_ha("left")
    ax.title.set_position([0.0,1.0])
    
    for i in range(2):
        ax = axs[i]
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.text(0, 1.15, chr(97 + i), transform=ax.transAxes, 
            fontsize=18, fontweight='bold', va='top', ha='right')
        ax.set_ylim(-0.25,1.25)

    plt.savefig(f'{figures_dir}/results_fig_4.png',dpi=300)
    plt.close()

def fig5(baseline):
    outcomes,_,base = get_data(baseline)

    base = base.groupby("Timestep").sum().reset_index()
    base['carbon_sequestered'] *= 1e-6 * (1/run_scale)
       
    plot = sb.lineplot(data=base,x="Timestep",y="carbon_sequestered")
    # assign_endlabels(plot,base,"Timestep","carbon_sequestered","Land use category",lpad=4)
    plot.set_xlabel("")
    plot.set_ylabel(r"Carbon sequestered (Mt CO$_{2}$/year)")
    # plt.subplots_adjust(right=0.8)
    plot.spines['right'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plot.set_ylim(0,80)
    plot.set_xlim(2025,2050)
    plot.axhline(y=35,linestyle=':',color='black')
    plot.text(2032, 41, 'Present day gross \n'+r'annual CO$_2$ emissions', horizontalalignment='center', verticalalignment='center')

    n2obase = 8
    ch4base = 30
    base['co2e_over_time'] = n2obase * outcomes['Nitrous oxide']/100 + ch4base * outcomes['Methane']/100 + 35
    sb.lineplot(data=base,x='Timestep',y='co2e_over_time',color='black',dashes=[1,1.5])
    plot.text(2032, 78, 'Below, plus modelled\n'+r'agricultural CO$_2$e emissions', horizontalalignment='center', verticalalignment='center')
    # plot.axhline(y=55,linestyle=':',color='black')
    # plot.text(2037, 59, 'Present day net annual emissions', horizontlalignment='center', verticalalignment='center')
    plt.savefig(f'{figures_dir}/results_fig_5.png',dpi=300)
    plt.close()
    return None # come back to carbon cost curve later
    # read in outcomes_aggregated from each small scale run, scale the carbon abatement by that in the $70 run, plot results
    carbondf = base[['Timestep','carbon_sequestered']].copy()
    carbondf['Initial carbon price'] = 70

    comprun = pd.read_csv(f"outputs/2024_03_22_18_00_00e0.0n0.0p0.0sed0.0n2o0.0ch40.0co270bfalsenolucfalsepathparabolic/plots/outcomes_aggregated.csv")
    scalar =  base.loc[base['Timestep']==2050,'carbon_sequestered'].sum() / (1e-6 * comprun.loc[comprun['Timestep']==25,'carbon_sequestered'].sum())

    for i in list(range(0,140,5)) + list(range(190,205,5)):
        run = pd.read_csv(f"outputs/2024_03_22_18_00_00e0.0n0.0p0.0sed0.0n2o0.0ch40.0co2{i}bfalsenolucfalsepathparabolic/plots/outcomes_aggregated.csv")
        print(i, run.loc[(run['Land use']=='carbon-forestry') & (run['Timestep']==25),'Area (ha)'].sum())
        run = run.groupby("Timestep").sum().reset_index()
        run['Timestep'] += 2025
        run['carbon_sequestered'] *= 1e-6 * scalar
        to_add = run[['Timestep','carbon_sequestered']].copy()
        to_add['Initial carbon price'] = i
        carbondf = pd.concat([carbondf,to_add])
    
    only_2050 = carbondf.loc[carbondf['Timestep']==2050]
    plot = sb.relplot(data=only_2050,y='Initial carbon price',x='carbon_sequestered',kind='scatter')
    for ax in plot.axes.flat:
        sb.regplot(y='Initial carbon price', x='carbon_sequestered', data=only_2050, ax=ax, scatter=False, lowess=True)
    plot.set_xlabels(r"2050 annual carbon sequestration (Mt CO$_{2}$/year)")
    plot.set_ylabels(r"Initial carbon price (NZD/t CO$_2$)")

    plt.savefig(f'{figures_dir}/results_fig_5b.png',dpi=300)
    

def rasterise_attribute(gdf,attribute,res=1000):
    bbox = gdf.total_bounds
    x_res = (bbox[2] - bbox[0]) / res
    y_res = (bbox[3] - bbox[1]) / res
    raster = rasterize([(geom,value) for geom, value in zip(gdf.geometry,gdf[attribute])], 
                    out_shape=(res,res), 
                    transform=rasterio.transform.from_origin(bbox[0], bbox[3], x_res, y_res),
                    fill=np.nan,
                    all_touched=True
                    )
    metadata = {
        'driver': 'GTiff',  # GeoTIFF format
        'dtype': 'int32',  # Data type of the raster; adjust as needed
        'nodata': None,  # Define the NoData value; adjust as needed
        'width': res,  # Number of columns of the raster
        'height': res,  # Number of rows of the raster
        'count': 1,  # Number of bands; typically 1 for a classified raster
        'crs': 'EPSG:2193',  # Coordinate Reference System
        'transform': rasterio.transform.from_origin(bbox[0], bbox[3], x_res, y_res)
    }

    # with rasterio.open(f'test.tif', 'w', **metadata) as dst:
    #     dst.write(raster, 1)
    return raster

def datafig(grid_cell):
    gdf = pyogrio.read_dataframe(grid_cell)
    
    fig, axs = plt.subplots(2,3,figsize=(10,8),gridspec_kw={'hspace':0.3})
    ax = axs[0,0]
    sat_img = plt.imread("grid_-37.7_174.8.png")
    ax.imshow(sat_img,extent=gdf.total_bounds,aspect="auto")
    ax.title.set_text("Satellite image")

    ax = axs[0,1]
    angle = rasterise_attribute(gdf,"angl_avg",res=1000)
    # palette = sb.choose_colorbrewer_palette("sequential",as_cmap=True)
    sb.heatmap(angle,vmin=0,ax=ax,cmap="YlOrRd")
    ax.title.set_text("Slope")

    ax = axs[0,1]

    ax = axs[0,2]
    gdf["DairyNZ revenue per ha"] = gdf["PastureYield_Rainfed200NCap_25m"] * 9454/16064
    dairy_rev = rasterise_attribute(gdf,"DairyNZ revenue per ha",res=1000)
    sb.heatmap(dairy_rev,ax=ax,cmap="YlGn")
    ax.title.set_text("Dairy revenue")

    ax = axs[1,0]
    max_n = rasterise_attribute(gdf,"dairy_typology_max_N",res=1000)
    sb.heatmap(max_n,ax=ax,cmap="BuPu")
    ax.title.set_text("Max N leaching")


    ax = axs[1,1]
    gdf["mayld_N"] = gdf["catchment_MAL_Nitrogen"] / (gdf["catchment_WaterShedAreaKM2"] * 100)
    mayld_n = rasterise_attribute(gdf,"mayld_N",res=1000)
    sb.heatmap(mayld_n,ax=ax,vmin=20,cmap="Reds")
    ax.title.set_text("Max allowed N yield")


    ax = axs[1,2]
    radiata = rasterise_attribute(gdf,"economic_indicator_revenue-cherry_value",res=1000)
    sb.heatmap(radiata,ax=ax,vmin=1e5,cmap="GnBu")
    ax.title.set_text("Cherry revenue")


    for i in range(6):
        ax = axs[i//3,i%3]
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.text(0, 1.2, chr(97 + i), transform=ax.transAxes, 
            fontsize=18, fontweight='bold', va='top', ha='right')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.title.set_ha("left")
        ax.title.set_position([0.0,1.0])
        # ax.set_ylim(-0.5,1.5)

    plt.savefig(f'{figures_dir}/data_fig.png',dpi=300)
    plt.close()


def process_farm_profits(land_use,dir):
    df = pd.read_feather(dir+'agents_0.arrow')
    by_agents = df.groupby(['uid','land_use_name']).sum().reset_index()
    subset_to_lu = by_agents[by_agents['land_use_name'] == land_use]
    subset_to_lu['profit_per_hectare'] = subset_to_lu['total_yearly_profit'] / subset_to_lu['area_ha']
    print(len(subset_to_lu), " farms in ", land_use)
    print(len(set(by_agents.loc[by_agents['land_use_name']!=land_use,'uid']) & set(subset_to_lu['uid'])), " agents managing other lus as well")
    return subset_to_lu['profit_per_hectare'].to_list()

def process_blnz_data():
    # datafile available from https://tools.beeflambnz.com/profitability-calculator
    blnz = pd.read_csv("data/B+LNZ export tool.csv")
    profit_intervals = np.arange(-300,1901,100)
    # trim the overflow bins
    freqs = blnz.iloc[0].to_list()
    rehistify = []
    for p, freq in zip(profit_intervals,freqs):
        for _ in range(int(freq*10)):
            rehistify.append(p+50)
    return rehistify

def snb_profit_fig(dir,suffix=""):
    land_use = "sheep-and-beef"
    model_data = process_farm_profits(land_use,dir)
    real_world_data= process_blnz_data()
    data = pd.DataFrame(zip(real_world_data,model_data),columns=["B+LNZ benchmark","MEIHAWAI"])
    melted = pd.melt(data)
    fig, ax = plt.subplots()
    # model_data = model_data[:100]
    
    sb.histplot(melted,x='value',ax=ax,stat="percent",bins=20,hue="variable")
    
    ax.set_xlabel("EBITr ($/ha/year)")
    ax.set_ylabel("Percentage of farms")
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    plt.savefig(f'{figures_dir}/profit_fig_{land_use}{suffix}.png',dpi=300)
    plt.close()


def process_dnz_eco_survey(count):
    pixel_height = {
        -250: 107,
        250: 144,
        750: 467,
        1250: 682,
        1750: 539,
        2250: 716,
        2750: 252,
        3250: 108,
        3750: 74,
        4250: 36,

    }
    histify = []
    tot = sum(pixel_height.values())
    for key,value in pixel_height.items():
        histify.extend([key]*int(count * value / tot))
    return histify
        

def dairy_profit_fig(dir,suffix=""):
    land_use = "dairy"
    model_data = process_farm_profits(land_use,dir)
    tax_return_data = np.random.normal(1909,1508,len(model_data))
    dnz_eco_survey_data = process_dnz_eco_survey(len(model_data))
    fig,axs = plt.subplots(2,1,figsize=(10,10),gridspec_kw={'hspace':0.3})
    dnz_bins = np.arange(-1500,6001,500)
    
    ax = axs[0]
    data = pd.DataFrame(zip(tax_return_data,model_data),columns=["Dairy tax returns","MEIHAWAI"])
    melted = pd.melt(data)
    sb.histplot(melted,x='value',ax=ax,stat="percent",bins=dnz_bins,hue="variable")
    ax.set_xlabel("")
    ax.set_ylabel("Percentage of farms")
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.text(0, 1.15, chr(97), transform=ax.transAxes, 
            fontsize=18, fontweight='bold', va='top', ha='right')

    ax = axs[1]
    data = pd.DataFrame(zip(dnz_eco_survey_data,model_data),columns=["DairyNZ economic survey","MEIHAWAI"])
    melted = pd.melt(data)
    sb.histplot(melted,x='value',ax=ax,stat="percent",bins=dnz_bins,hue="variable")
    ax.set_xlabel("EBITr ($/ha/year)")
    ax.set_ylabel("Percentage of farms")
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.text(0, 1.15, chr(98), transform=ax.transAxes, 
        fontsize=18, fontweight='bold', va='top', ha='right')

        
    plt.savefig(f'{figures_dir}/profit_fig_{land_use}_{suffix}.png',dpi=300)
    plt.close()

def dairy_snb_carbon_profit_fig(dir):
    dairy_prof = process_farm_profits("dairy",dir)
    snb_prof = process_farm_profits("sheep-and-beef",dir)
    fig,ax = plt.subplots(gridspec_kw={'bottom':0.15})
    data = pd.DataFrame(zip(snb_prof,dairy_prof),columns=["Sheep and beef","Dairy"])
    melted = pd.melt(data)
    plot = sb.histplot(melted,x='value',ax=ax,stat="percent",bins=20,hue="variable",legend=True)
    ax.set_xlabel("Operating profit ($/ha/year)")
    ax.set_ylabel("Percentage of farms")
    # unbelievable
    plot.legend(['Dairy','Sheep and beef'],reverse=True,title="",bbox_to_anchor=(0.5,1.1), borderaxespad=0.1,loc='upper center',ncol=2)
    carbon_npv = 328.62
    annuity_scalar = 1/15.81
    carbon20 = carbon_npv*annuity_scalar*20
    ax.axvline(x=carbon20,ymax=0.8,linestyle=':',linewidth=2,color='black')
    ax.text(carbon20-700, 13, 'Carbon forestry,\nANPV ($20/ton)', horizontalalignment='center', verticalalignment='center',fontdict={'size':12})
    carbon70 = carbon_npv*annuity_scalar*70
    ax.axvline(x=carbon70,ymax=0.8,linestyle=':',linewidth=2,color='black')
    ax.text(carbon70, 14, 'Carbon forestry,\nANPV ($70/ton)', horizontalalignment='center', verticalalignment='center',fontdict={'size':12})
    carbon100 = carbon_npv*annuity_scalar*100
    ax.axvline(x=carbon100,ymax=0.8,linestyle=':',linewidth=2,color='black')
    ax.text(carbon100+800, 12, 'Carbon forestry,\nANPV ($100/ton)', horizontalalignment='center', verticalalignment='center',fontdict={'size':12})
    ax.set_ylim(0,16)
    

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    plt.savefig(f'{figures_dir}/profit_fig_dairy_snb.png',dpi=300)
    plt.close()


def fig6(dir):
    catchment_data = pd.read_feather(dir.replace('plots/','')+'catchments_0.arrow')
    
    catchment_data['N_excess'] = catchment_data['Load_N'] - catchment_data['MAL_N']
    positive_N = catchment_data[(catchment_data['N_excess']>0) & (catchment_data['MAL_N']>0)]
    positive_N['total_price'] = positive_N['N_excess'] * 100
    positive_N['price_per_kg'] = positive_N['total_price'] / positive_N['Load_N']
    print(positive_N['price_per_kg'].mean())
    plot = sb.displot(positive_N,x='price_per_kg',kind='hist',bins=range(0,101,10))
    plot.set_xlabels("Price per kg N ($)")
    plot.set_ylabels("Number of catchments")
    plt.savefig(f'{figures_dir}/results_fig_6a.png',dpi=300)
    plt.close()
    catchment_data['P_excess'] = catchment_data['Load_P'] - catchment_data['MAL_P']

    
    positive_P = catchment_data[(catchment_data['P_excess']>0) & (catchment_data['MAL_P']>0)]
    positive_P['total_price'] = positive_P['P_excess'] * 1500
    positive_P['price_per_kg'] = positive_P['total_price'] / positive_P['Load_P']
    print(positive_P['price_per_kg'].mean())

    plot = sb.displot(positive_P,x='price_per_kg',kind='hist',bins=range(0,1501,100))
    plot.set_xlabels("Price per kg P ($)")
    plot.set_ylabels("Number of catchments")
    plt.savefig(f'{figures_dir}/results_fig_6b.png',dpi=300)
    plt.close()


def progressive_pricing_fig():
    df = pd.read_csv('outputs/nesiparamsearch.csv')
    df = df[df['include']==True]
    melted = pd.melt(df,id_vars=['25-year profit deficit (rel)'],value_vars=['N target','P target','CH4 target','N2O target','Sediment target'])
    rename_vars = {
        'N target':'Nitrogen',
        'P target':'Phosphorous',
        'CH4 target':'Methane',
        'N2O target':'Nitrous oxide',
        'Sediment target':'Sediment'
    }
    melted['variable'] = melted['variable'].apply(lambda x: rename_vars[x])
    melted['value'] *= 100
    melted['25-year profit deficit (rel)'] *= 100
    melted['Target met'] = melted.apply(lambda x: x['value'] < 65 if x['variable'] in ['Methane','Nitrous oxide'] else x['value'] < 10,axis=1)
    melted.rename(columns={'variable':'Pollutant'},inplace=True)
    fig, ax = plt.subplots(1,1,figsize=(8,5),gridspec_kw={'top':0.8})
    sb.lineplot(data=melted,x='25-year profit deficit (rel)',y='value',hue='Pollutant',ax=ax,size=1,legend=False)#,size=1)
    sb.scatterplot(data=melted,x='25-year profit deficit (rel)',y='value',hue='Pollutant',ax=ax,style='Target met',s=50)
    ax.set_xlabel("25-year decrease in producer profit (%)")
    ax.set_ylabel("2050 pollution level index (2025 = 100)")
    ax.set_ylim(0,100)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.axhline(y=65,linestyle=':',color='black')
    ax.text(9, 60, r'CH$_4$, N$_2$O target', fontname='Arial', horizontalalignment='center', verticalalignment='center')
    ax.axhline(y=10,linestyle=':',color='black')
    ax.text(9, 15, 'N, P, Sediment target', horizontalalignment='center', verticalalignment='center')
    ax.legend(title="",bbox_to_anchor=(0.5,1.2), frameon=False, borderpad=0.2,loc='upper center',ncol=3,fontsize=14,labelspacing=0.2)
    plt.savefig(f'{figures_dir}/results_fig_7.png',dpi=300)

    plt.close()



fig1(baseline_dir,fully_optimised_dir)
fig2v2(baseline_dir)
fig3(baseline_dir,fully_optimised_dir)
fig3(baseline_dir,zerocarbon,suffix="_zerocarbon")
fig3(baseline_dir,extensiononlylow,suffix="_extensiononlylow")
fig4(baseline_dir,fully_optimised_dir)
fig4aares(fully_optimised_dir)
fig5(baseline_dir)

