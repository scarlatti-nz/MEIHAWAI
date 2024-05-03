import matplotlib
import seaborn as sb
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scipy.stats as stats
import random
import os
import argparse
import numpy as np
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import importlib
import math
# sb.set_font('Arial')
sb.set_palette("bright")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 14
A=5300
B=1000
C=2000
D=74 # max nitrogen loss
E=1000 # mitigation cost
alpha=26.9
b=73.75
c=73.75
gamma=0.52

CH4 = 8000
N2O = 2826

def prod_int(x,cap,h=20):
    return A*(1-np.e**(-x * h * cap))

def prof_int(x,cap,vaccine=False,h=20):
    cost_adjuster = 300 if vaccine else 0
    return A*(1-np.e**(-x * h * cap)) - B * x - C - cost_adjuster

def nl_mit(x,mcap,intensity=0.3,h=10):
    # intensity = 0.5
    return  (1 - gamma*(1-np.e**(-x*mcap*h))) * (alpha + b * intensity + c * intensity**2)
def nl_int(x, mcap, mitigation=0.3,h=10):
    # mitigation = 0.5
    min_nl = (alpha + b * 0.24 + c * 0.24**2)
    max_nl = (alpha + b * 0.4 + c * 0.4**2)
    actual = (alpha + b * x + c * x**2)
    premit = np.clip(actual,min_nl,max_nl)
    return (1 -  gamma*(1-np.e**(-mitigation * mcap*h))) * premit

def prof_mit(x,intensity=0.5,cap=0.3,h=20):
    return A*(1-np.e**(-intensity * h * cap)) - B*intensity - C - E * x

def ghg_int(x,cap,vaccine=False,h=20):
    methane = CH4 * (1 - np.e**(-x * h * cap)) * (0.3 if vaccine else 1)
    nitrous_oxide = N2O * (alpha + b * x + (1 - alpha - b) * x**2)
    return methane + nitrous_oxide

def figs_for_docs():
    # Define the x values for your four plots
    x = np.linspace(0, 1, 100)
    eqs = [prod_int, prof_int, nl_mit, nl_int]
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    pcaps = [0.4,0.55,0.7]
    mcaps = [0.08,0.4,0.72]
    caplabels = ["10th percentile", "50th percentile", "90th percentile"]
    yaxes = ["Revenue ($/ha/yr)","Profit ($/ha/yr)","Nitrogen loss (kg/ha/yr)","Nitrogen loss (kg/ha/yr)"]
    xaxes = ["Intensity (I)","Intensity (I)","Mitigation (M)","Intensity (I)"]
    ylims = [(4000,5600),(1400,3000),(20,100),(20,100)]
    xlims = [(0.2,0.5),(0.2,0.5),(0,1),(0.2,0.5)]
    tickdist = [0.05,0.05,0.2,0.05]
    anspace = [(),(0.3,0.45),(0.7,0.55),(0.7,0.55)]
    annotations = ["","Mitigation = 0","Intensity = 0.3","Mitigation = 0.3"]
    # Plot the three series for each subplot
    for i in range(4):
        ltitle = "Production capability" if i < 2 else "Mitigation capability"
        if i<2:
            caps = pcaps
        else:
            caps = mcaps
        eq = eqs[i]
        y1 = eq(x,caps[0])
        y2 = eq(x,caps[1])
        y3 = eq(x,caps[2])
        ax = axs[i//2, i%2]
        ax.plot(x, y1, label=str(caplabels[0]))
        ax.plot(x, y2, label=str(caplabels[1]))
        ax.plot(x, y3, label=str(caplabels[2]))
        ax.set(xlabel=xaxes[i], ylabel=yaxes[i], ylim=ylims[i], xlim=xlims[i]) 
        ax.text(-0.1, 1.1, chr(97 + i), transform=ax.transAxes, 
            fontsize=14, fontweight='bold', va='top', ha='right')
        ax.legend(title=ltitle)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_major_locator(MultipleLocator(tickdist[i]))
        ax.tick_params(which='both', direction='out', length=6)
        # add annotations
        annotation_text = annotations[i]
        if annotation_text != "":
            ax.text(anspace[i][0], anspace[i][1], annotation_text, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    plt.savefig(f'plots/figs_for_docs_2.png',dpi=300)
    plt.close()

    ## outcome curves

    ## LBD
    # def lbd(x0,Q,lr,ts):
    #     x = [x0]
    #     for i in range(ts):
    #         x.append(x[i] + Q*lr*(1-x[i]))
    #     return x
    # fig,ax = plt.subplots(1,2,figsize=(10,10))
    # x = np.arange(0,41)
    # y1 = lbd(0.1,0.08,0.6,40)
    # y2 = lbd(0.1,0.08,0.8,40)
    # y3 = lbd(0.1,0.08,1,40)

    # ax[0].plot(x,y1,label="Learning rate = 0.6")
    # ax[0].plot(x,y2,label="Learning rate = 0.8")
    # ax[0].plot(x,y3,label="Learning rate = 1")
    # B = 3500
    # C = 1000
    # def prof_cap_optimal_intensity(cap):
    #     x = max(0,min(1,math.log((A*20*cap/B),np.e)/(20*cap)))
    #     print(cap,x)
    #     return prof_int(x,cap)
    # z1 = [prof_cap_optimal_intensity(y) for y in y1]
    # z2 = [prof_cap_optimal_intensity(y) for y in y2]
    # z3 = [prof_cap_optimal_intensity(y) for y in y3]
    # ax[1].plot(x,z1,label="Learning rate = 0.6")
    # ax[1].plot(x,z2,label="Learning rate = 0.8")
    # ax[1].plot(x,z3,label="Learning rate = 1")
    # plt.show()
figs_for_docs()

def outcome_curves():
    I = np.linspace(0,1,100)
    M = np.linspace(0,1,100)
    
    fig, axs = plt.subplots(2,2,figsize=(10,10))
    axs[0,0].plot(nl_int(I,0.5,mitigation=0),prof_int(I,0.8),label="0.8")
    axs[0,0].plot(nl_int(I,0.5,mitigation=0),prof_int(I,0.4),label="0.4")
    axs[0,0].plot(nl_int(I,0.5,mitigation=0),prof_int(I,0.2),label="0.2")
    
    axs[0,0].set(xlabel="Nitrogen loss (kg/ha/yr)",ylabel="Profit ($/ha/yr)")
    axs[0,0].legend(title="Production capability")
    axs[0,0].annotate('',xy=(60, 2700), xytext=(90, 2450), arrowprops=dict(color='black', arrowstyle='<->'))
    axs[0,0].text(80, 2800, 'Varying intensity', horizontalalignment='center', verticalalignment='center')
    
    # set pc=0.3 and 
    axs[0,1].plot(nl_mit(M,0.8,intensity=0.7),prof_mit(M,intensity=0.7,cap=0.3),label="0.8")
    axs[0,1].plot(nl_mit(M,0.4,intensity=0.7),prof_mit(M,intensity=0.7,cap=0.3),label="0.4")
    axs[0,1].plot(nl_mit(M,0.2,intensity=0.7),prof_mit(M,intensity=0.7,cap=0.3),label="0.2")
    
    axs[0,1].set(xlabel="Nitrogen loss (kg/ha/yr)",ylabel="Profit ($/ha/yr)")
    axs[0,1].legend(title="Mitigation capability")
    axs[0,1].annotate('',xy=(50, 2550), xytext=(60, 2600), arrowprops=dict(color='black', arrowstyle='<->'))
    axs[0,1].text(55, 2750, 'Varying mitigation', horizontalalignment='center', verticalalignment='center')

    #ghgs
    axs[1,0].plot(ghg_int(I,0.8),prof_int(I,0.8),label="0.8")
    axs[1,0].plot(ghg_int(I,0.4),prof_int(I,0.4),label="0.4")
    axs[1,0].plot(ghg_int(I,0.2),prof_int(I,0.2),label="0.2")
    
    axs[1,0].set(xlabel="GHG emissions (kg CO2e/ha/yr)",ylabel="Profit ($/ha/yr)")
    axs[1,0].legend(title="Production capability")
    axs[1,0].annotate('',xy=(5300,1000), xytext=(7100, 2000), arrowprops=dict(color='black', arrowstyle='<->'))
    axs[1,0].text(8560, 1500, 'Varying intensity', horizontalalignment='center', verticalalignment='center')

    #ghg pricing
    axs[1,1].plot(ghg_int(I,0.4,vaccine=True),prof_int(I,0.4,vaccine=True),label="Used")
    axs[1,1].plot(ghg_int(I,0.4),prof_int(I,0.4),label="Not used")
    axs[1,1].set(xlabel="GHG emissions (kg CO2e/ha/yr)",ylabel="Profit ($/ha/yr)")
    axs[1,1].legend(title="Bromoform bolus")
    axs[1,1].annotate('',xy=(5000,1000), xytext=(6500, 2000), arrowprops=dict(color='black', arrowstyle='<->'))
    axs[1,1].text(8000, 1500, 'Varying intensity', horizontalalignment='center', verticalalignment='center')

    for i in range(4):
        ax = axs[i//2, i%2]
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.text(-0.1, 1.1, chr(65 + i), transform=ax.transAxes, 
            fontsize=14, fontweight='bold', va='top', ha='right')
        ax.set_ylim(0,3000)
        # reverse axis legend order
        
    plt.subplots_adjust(wspace=0.4, hspace=0.3)
    plt.savefig(f'plots/outcome_curves_1.png',dpi=300)
    plt.close()
outcome_curves()



def prod_cap_error():
    data = pd.read_csv('outputs/production capability error.csv')
    fig,axes = plt.subplots(2,1,figsize=(10,10),sharex=True,gridspec_kw={'height_ratios': [5, 1]})
    ax = axes[0]
    ## main plot
    
    x = data[data['pc']==0.3]['pce']
    y1 = data[data['pc']==0.2]['util']
    y2 = data[data['pc']==0.3]['util']
    y3 = data[data['pc']==0.5]['util']
    ax.plot(x,y1,label="0.2")
    ax.plot(x,y2,label="0.3")
    ax.plot(x,y3,label="0.5")
    ax.set(ylabel="True utility ANPV ($/ha/yr)")
    ax.legend(title="Production capability")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xscale('log')
    ax.set_ylim(0,500)
    ax = axes[1]
    ## error plot
    x = data[data['env_scalar']==10]['pce']
    y1 = data[data['env_scalar']==10]['intensity']
    ax.plot(x,y1,label="Very high")
    ax.set(xlabel="Deviation parameter (log scale)")
    x = np.linspace(0.3,3,100)
    def density(x,mu,sigma):
        return 1/(x*sigma*np.sqrt(2*np.pi))*np.exp(-((np.log(x)-mu)**2)/(2*sigma**2))
    y_density = density(x,-0.1,0.7)
    ax.plot(x,y_density)
    ax.set(ylabel="Density")
    ax.axvline(x=np.e**(-0.1),color="black",ymin=0,ymax=0.72,linestyle="-")
    ax.annotate('Median', xy=(np.e**(-0.1), 0.5), xytext=(0.7, 0.4),
             arrowprops=dict(facecolor='black', arrowstyle='-'))
    plt.yticks([],[])
    plt.xticks([0.3,0.5,0.7,1,1.5,2,3],["0.3","0.5","0.7","1","1.5","2","3"])
    plt.savefig(f'plots/prod_cap_error.png',dpi=300)
    plt.close()
prod_cap_error()