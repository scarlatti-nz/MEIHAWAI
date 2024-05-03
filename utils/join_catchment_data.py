import pandas as pd
import pyogrio
import matplotlib
import sys

matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['font.size'] = 14

import matplotlib.pyplot as plt
import scipy as sp
from scipy import stats
import geopandas as gpd
import seaborn as sb
sb.set_palette("bright")

# catchment_datafile = '../nesi_downloads/meihawai_outputs/2024_03_23_00_28_42e0.0n0.0p0.0sed0.0n2o0.0ch40.0co270.0bfalsenolucfalsepathparabolic/plots/catchment_data.csv'
# catchment_data = pd.read_csv(catchment_datafile)

catchment_datafile_sediment_default = '../nesi_downloads/meihawai_outputs/2024_03_27_15_21_08e0.0n0.0p0.0sed0.0n2o0.0ch40.0co270.0bfalsenolucfalsepathparabolic/catchments_1.arrow'
# catchment_datafile_sediment_alt = '../nesi_downloads/meihawai_outputs/2024_04_04_18_21_02e0.0n0.0p0.0sed0.0n2o0.0ch40.0co270.0bfalsenolucfalsepathparabolic/catchments_1.arrow'
catchment_datafile_sediment_alt = '../nesi_downloads/catchments_0.arrow'

if 'alt' in sys.argv:
    catchment_data = pd.read_feather(catchment_datafile_sediment_alt)
else:
    catchment_data = pd.read_feather(catchment_datafile_sediment_default)

catchment_data = catchment_data.rename(columns={'nzseg':'nzsegment'})


exempt_segs = pd.unique(catchment_data[catchment_data['has_Sediment_exemption'] == True]['nzsegment'])

print(len(exempt_segs)/len(pd.unique(catchment_data['nzsegment'])) * 100, '% of catchments have a sediment exemption')
# catchment_data = catchment_data[catchment_data['Timestep'] == 0]
catchment_data = catchment_data[catchment_data['nzsegment'] != 0]
# catchment_data = catchment_data[~ catchment_data['nzsegment'].isin(exempt_segs)]

# data files can be downloaded from whitiwhiti ora:
# shapefile - https://landuseopportunities.nz/dataset/current-state-vs-freshwater-objectives-nutrients-nitrogen-and-phosphorus/resource/f5beac74-4b65-4be4-a870-5236def92609
# MAL csvs - not yet on whitiwhiti ora, available on request.
shapefile = pyogrio.read_dataframe("data/watershed_nutrients_v1-aug-2023-MAL.shp")
MAL_csv = pd.read_csv("data/MAL and catchment data.csv")
sediment_MAL_csv = pd.read_csv("data/Sediment25Nov2023_RECout.csv")

joined = shapefile.merge(catchment_data, on='nzsegment', how='right')

joined2 = joined.merge(MAL_csv, on='nzsegment', how='left')
joined2 = joined2.merge(sediment_MAL_csv, on='nzsegment', how='left', suffixes=('', '_Sediment'))
joined2 = joined2[joined2['total_coverage'] < 1.5]

# joined2['Load_Sediment'] = joined2['Load_Sediment'] * 315


joined2['Yld_N'] = joined2['Load_N'] / (100 * joined2['WaterShedAreaKM2'])
joined2['ExYld_N'] = (joined2['Load_N'] - joined2['MAL_N']) / (100 * joined2['WaterShedAreaKM2'])
joined2['Yld_P'] = joined2['Load_P'] / (100 * joined2['WaterShedAreaKM2'])
joined2['ExYld_P'] = (joined2['Load_P'] - joined2['MAL_P']) / (100 * joined2['WaterShedAreaKM2'])
joined2['Yld_Sediment'] = joined2['Load_Sediment'] / (100 * joined2['WaterShedAreaKM2'])
joined2['ExYld_Sediment'] = (joined2['Load_Sediment'] - joined2['MAL_Sediment']) / (100 * joined2['WaterShedAreaKM2'])
joined2['CurrentSedYield'] = joined2['CurrentSedYield'] / 100
joined2['MAYLD_Sediment'] = joined2['MAL_Sediment'] / (100 * joined2['WaterShedAreaKM2'])


# report figures here
figures_dir = "figures/"

N_to_hist = joined2[['TNYield','Yld_N']]
N_to_hist = N_to_hist.melt()
N_to_hist = N_to_hist[N_to_hist['value'] < 40]
plot = sb.displot(N_to_hist, x='value', hue='variable', kind='hist',bins=15,legend=False)
plt.annotate('MEIHAWAI', xy=(5.3,80000), xytext=(15,90000), 
            xycoords='data', textcoords = 'data', 
            ha='center', va='center', fontsize=14,
            arrowprops = {
                'arrowstyle': '-',
                'color': 'black'
            })
plt.annotate('Snelder et al.', xy=(8,55000), xytext=(18,65000),
            xycoords='data', textcoords = 'data',
            ha='center', va='center', fontsize=14,
            arrowprops = {
                'arrowstyle': '-',
                'color': 'black'
            })
plt.xlabel('Area-normalised N load (kg/ha/year)')
plt.ylabel('DN2.4 segments')
plt.savefig(figures_dir + '/N_yield_hist.png',dpi=300)
plt.close()

P_to_hist = joined2[['TPYield','Yld_P']]
P_to_hist = P_to_hist.melt()
P_to_hist = P_to_hist[P_to_hist['value'] < 2]
plot = sb.displot(P_to_hist, x='value', hue='variable', kind='hist',bins=15,legend=False)
plt.annotate('MEIHAWAI', xy=(0.13,80000), xytext=(0.9,90000),
            xycoords='data', textcoords = 'data',
            ha='center', va='center', fontsize=14,
            arrowprops = {
                'arrowstyle': '-',
                'color': 'black'
            })
plt.annotate('Snelder et al.', xy=(0.4,100000), xytext=(1.1,110000),
            xycoords='data', textcoords = 'data',
            ha='center', va='center', fontsize=14,
            arrowprops = {
                'arrowstyle': '-',
                'color': 'black'
            })
plt.xlabel('Area-normalised P load (kg/ha/year)')
plt.ylabel('DN2.4 segments')
plt.savefig(figures_dir + '/P_yield_hist.png',dpi=300)
plt.close()

Sed_to_hist = joined2[['CurrentSedYield','Yld_Sediment']]
Sed_to_hist = Sed_to_hist.melt()

quintiles = Sed_to_hist.loc[Sed_to_hist['variable'] == 'Yld_Sediment','value'].quantile([0.2,0.4,0.6,0.8,1.0])
quintiles = quintiles.values.tolist()
bardata = pd.DataFrame(columns=['bin','model','count'])

fig,ax = plt.subplots(2,1,figsize=(10,6),gridspec_kw={'hspace':0.3})
for l,q in zip([0]+quintiles[:4],quintiles):
    
    bardata.loc[len(bardata)] = [f'{l:.2f} - {q:.2f}', 'Snelder et al.', len(Sed_to_hist[(Sed_to_hist['value'] <= q) & (Sed_to_hist['value']>l) & (Sed_to_hist['variable'] == 'CurrentSedYield')])]
    bardata.loc[len(bardata)] = [f'{l:.2f} - {q:.2f}', 'MEIHAWAI', len(Sed_to_hist[(Sed_to_hist['value'] <= q) & (Sed_to_hist['value']>l) & (Sed_to_hist['variable'] == 'Yld_Sediment')])]
sb.barplot(data=bardata, x='bin', y='count', hue='model', ax=ax[0])
ax[0].set_xlabel('')
ax[0].set_ylabel('DN2.4 segments')
ax[0].legend(title='', loc='upper center', bbox_to_anchor=(0.5, 1.3), ncol=2, title_fontsize='large')
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[0].text(0, 1.15, chr(97), transform=ax[0].transAxes, 
            fontsize=18, fontweight='bold', va='top', ha='right')
# upperquintiles = Sed_to_hist.loc[(Sed_to_hist['value'] > quintiles[3]) & (Sed_to_hist['variable'] == 'Yld_Sediment'),'value'].quantile([0.2,0.4,0.6,0.8,1.0])
# upperquintiles = upperquintiles.values.tolist()
upperquintiles = [5,10,25,100,quintiles[4]]

bardata2 = pd.DataFrame(columns=['bin','model','count'])
for l,q in zip([quintiles[3]]+upperquintiles[:4],upperquintiles):
    bardata2.loc[len(bardata2)] = [f'{l:.2f} - {q:.2f}', 'Snelder et al.', len(Sed_to_hist[(Sed_to_hist['value'] <= q) & (Sed_to_hist['value']>l) & (Sed_to_hist['variable'] == 'CurrentSedYield')])]
    bardata2.loc[len(bardata2)] = [f'{l:.2f} - {q:.2f}', 'MEIHAWAI', len(Sed_to_hist[(Sed_to_hist['value'] <= q) & (Sed_to_hist['value']>l) & (Sed_to_hist['variable'] == 'Yld_Sediment')])]
sb.barplot(data=bardata2, x='bin', y='count', hue='model', ax=ax[1])


ax[1].set_xlabel('Area-normalised sediment load (t/ha/year)')
ax[1].set_ylabel('DN2.4 segments')
ax[1].get_legend().remove()
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].text(0, 1.15, chr(98), transform=ax[1].transAxes, 
            fontsize=18, fontweight='bold', va='top', ha='right')

if 'alt' in sys.argv:
    plt.savefig(figures_dir + '/Sediment_yield_barplot_alt.png',dpi=300)
else:
    plt.savefig(figures_dir + '/Sediment_yield_barplot.png',dpi=300)
plt.close()


# Sed_to_hist = Sed_to_hist[Sed_to_hist['value'] < 10]
# Sed_to_hist = Sed_to_hist[Sed_to_hist['value'] > 0.1]
# plot = sb.displot(Sed_to_hist, x='value', hue='variable', kind='hist',bins=15,legend=False,log_scale=True)
# plt.yscale('log')
# plt.annotate('MEIHAWAI', xy=(0.67,140000), xytext=(3,120000),
#             xycoords='data', textcoords = 'data',
#             ha='center', va='center', fontsize=14,
#             arrowprops = {
#                 'arrowstyle': '-',
#                 'color': 'black'
#             })
# plt.annotate('Snelder et al.', xy=(1.32,54000), xytext=(3.5,65000),
#             xycoords='data', textcoords = 'data',
#             ha='center', va='center', fontsize=14,
#             arrowprops = {
#                 'arrowstyle': '-',
#                 'color': 'black'
#             })
# plt.xlabel('Area-normalised sediment load (t/ha/year)')
# plt.ylabel('DN2.4 segments')
# plt.savefig(figures_dir + '/Sediment_yield_hist.png',dpi=300)
# plt.close()

# pyogrio.write_dataframe(joined2, catchment_datafile.replace('.csv','.gpkg'))
joined2 = joined2[joined2['MAL_N'] > 0]
proportion_zero_N = len(joined2[joined2['excess_N_percent'] == 0]) / len(joined2)
joined2 = joined2[joined2['MAL_P'] > 0]
proportion_zero_P = len(joined2[joined2['excess_P_percent'] == 0]) / len(joined2)
proportion_zero_Sediment =  len(joined2[joined2['excess_Sediment_percent'] == 0]) / len(joined2)



# proportion_zero_P_ton = len(joined[joined['TP_ExYld%'] == 0]) / len(joined)
# proportion_zero_ton = len(joined[joined['TN_ExYld%'] == 0]) / len(joined)

joined2['ton_marginal_excess_N'] = joined2['TNYield'] - (joined2['MAL_N'] / joined2['watershed_area_ha'])
joined2['ton_marginal_excess_P'] = joined2['TPYield'] - (joined2['MAL_P'] / joined2['watershed_area_ha'])
joined2['ton_marginal_excess_Sediment'] = joined2['CurrentSedYield'] - (joined2['MAL_Sediment'] / joined2['watershed_area_ha'])

proportion_zero_N_ton = len(joined2[joined2['ton_marginal_excess_N'] < 0] ) / len(joined2)
proportion_zero_P_ton = len(joined2[joined2['ton_marginal_excess_P'] < 0] ) / len(joined2)
proportion_zero_Sediment_ton = len(joined2[joined2['ton_marginal_excess_Sediment'] < 0] ) / len(joined2)

print("Our model predicts that {}% of catchments will have zero marginal excess N in 2023".format(proportion_zero_N*100))
print("Ton's model predicts that {}% of catchments will have zero marginal excess N in 2023".format(proportion_zero_N_ton*100))

print('\n') 

print("Our model predicts that {}% of catchments will have zero marginal excess P in 2023".format(proportion_zero_P*100))
print("Ton's model predicts that {}% of catchments will have zero marginal excess P in 2023".format(proportion_zero_P_ton*100))

print('\n')

print("Our model predicts that {}% of catchments will have zero marginal excess sediment in 2023".format(proportion_zero_Sediment*100))
print("Ton's model predicts that {}% of catchments will have zero marginal excess sediment in 2023".format(proportion_zero_Sediment_ton*100))

print('\n')

# correl = stats.pearsonr(joined2['ton_marginal_excess_N'], joined2['ExYld_N'])
# print("The correlation between our model and Ton's model's marginal excess N yields is {}".format(correl[0]))

plt.scatter(joined2['ton_marginal_excess_N'], joined2['ExYld_N'], s=1)
plt.xlim(-50,50)
plt.ylim(-50,50)
plt.savefig('plots/excess_N_scatter.png')
plt.close()

plt.scatter(joined2['ton_marginal_excess_P'], joined2['ExYld_P'], s=1)
plt.xlim(-5,5)
plt.ylim(-5,5)
plt.savefig('plots/excess_P_scatter.png')
plt.close()

plt.scatter(joined2['ton_marginal_excess_Sediment'], joined2['ExYld_Sediment'], s=1)
plt.xlim(-25,25)
plt.ylim(-25,25)
plt.savefig('plots/excess_Sediment_scatter.png')
plt.close()


plt.scatter(joined2['TNYield'], joined2['Yld_N'], s=1)
plt.savefig('plots/N_yld_scatter.png')
plt.close()

plt.scatter(joined2['TPYield'], joined2['Yld_P'], s=1)
plt.savefig('plots/P_yld_scatter.png')
plt.close()

plt.scatter(joined2['CurrentSedYield'], joined2['Yld_Sediment'], s=1)
plt.savefig('plots/Sediment_yld_scatter.png')
plt.close()

# nonzero_us = joined[joined['excess_N_percent'] != 0]

# nonzero_ton = joined[joined['TN_ExYld%'] != 0]




# ton_rounded = nonzero_ton['TN_ExYld%'].apply(lambda x: round(x, 15))
# top_vals = ton_rounded.value_counts().head(5)
# print(top_vals)
# print(len(nonzero_ton))
# nonzero_ton = nonzero_ton[~ nonzero_ton['TN_ExYld%'].isin(top_vals.index)]
# print(len(nonzero_ton))


# plt.hist(nonzero_us['excess_N_percent'], bins=10)
# plt.savefig('plots/model_excess_N_hist.png')
# plt.close()

# plt.hist(nonzero_ton['TN_ExYld%'], bins=10)
# plt.savefig('plots/ton_excess_N_hist.png')
# plt.close()

# joined2 = joined2[joined2['Load_N'] / (100 * joined2['WaterShedAreaKM2']) < 20]

print(f"Ton's model predicts a mean yield of {joined2['TPYield'].mean()} kg P/ha. Standard deviation is {joined2['TPYield'].std()}")
print(f"Our model predicts a mean yield of {joined2['Yld_P'].mean()} kg P/ha. Standard deviation is {joined2['Yld_P'].std()}")

print('\n')

print(f"Ton's model predicts a mean yield of {joined2['TNYield'].mean()} kg N/ha. Standard deviation is {joined2['TNYield'].std()}")
print(f"Our model predicts a mean yield of {joined2['Yld_N'].mean()} kg N/ha. Standard deviation is {joined2['Yld_N'].std()}")

print('\n')

print(f"Ton's model predicts a mean yield of {joined2['CurrentSedYield'].mean()} t/ha. Median is {joined2['CurrentSedYield'].median()} t/ha. Standard deviation is {joined2['CurrentSedYield'].std()}")
print(f"Our model predicts a mean yield of {joined2['Yld_Sediment'].mean()} t/ha. Median is {joined2['Yld_Sediment'].median()} t/ha. Standard deviation is {joined2['Yld_Sediment'].std()}")

print('\n')



print(f"Ton's model predicts a mean marginal excess yield of {joined2['ton_marginal_excess_N'].mean()} kg N/ha. Standard deviation is {joined2['ton_marginal_excess_N'].std()}")
print(f"Our model predicts a mean marginal excess yield of {joined2['ExYld_N'].mean()} kg N/ha. Standard deviation is {joined2['ExYld_N'].std()}")

print('\n')

print(f"Ton's model predicts a mean marginal excess yield of {joined2['ton_marginal_excess_P'].mean()} kg P/ha. Standard deviation is {joined2['ton_marginal_excess_P'].std()}")
print(f"Our model predicts a mean marginal excess yield of {joined2['ExYld_P'].mean()} kg P/ha. Standard deviation is {joined2['ExYld_P'].std()}")

print('\n')

print(f"Ton's model predicts a mean marginal excess yield of {joined2['ton_marginal_excess_Sediment'].mean()} t/ha. Standard deviation is {joined2['ton_marginal_excess_Sediment'].std()}")
print(f"Our model predicts a mean marginal excess yield of {joined2['ExYld_Sediment'].mean()} t/ha. Standard deviation is {joined2['ExYld_Sediment'].std()}")







# plt.hist(joined2['Load_N'] / (100 * joined2['WaterShedAreaKM2']), bins=50)
# plt.savefig('plots/N_yld_hist.png')
# plt.close()

# plt.hist(joined2['ExYld_N'], bins=50)
# plt.savefig('plots/excess_N_hist.png')
# plt.close()

# plt.hist(joined2['ExYld_P'], bins=50)
# plt.savefig('plots/excess_P_hist.png')
# plt.close()

# plt.hist(joined2['ton_marginal_excess_P'], bins=50)
# plt.savefig('plots/ton_excess_P_hist.png')
# plt.close()

# plt.hist(joined2['ton_marginal_excess_N'], bins=50)
# plt.savefig('plots/ton_excess_N_hist.png')
# plt.close()


# plt.hist(joined2['TNYield'], bins=50)
# plt.savefig('plots/ton_N_yld_hist.png')
# plt.close()

# plt.hist(joined2['MAL_N'] / (100 * joined2['WaterShedAreaKM2']), bins=20)
# plt.savefig('plots/our_MAYld_hist.png')
# plt.close()


# plt.hist(joined2['Yld_Sediment'], bins=50)
# plt.savefig('plots/Sediment_yld_hist.png')
# plt.close()

# plt.hist(joined2['CurrentSedYield'], bins=50)
# plt.savefig('plots/ton_Sediment_yld_hist.png')
# plt.close()



# plt.hist(joined2.loc[joined2['isTerminal']==True,'ton_MAyld_method_2'], bins=20)
# plt.savefig('plots/ton_MAYld_hist.png')
# plt.close()

# plt.scatter(joined2.loc[joined2['isTerminal']==True,'ton_MAyld_method_2'], joined2.loc[joined2['isTerminal']==True,'MAL_N'] / (100 * joined2.loc[joined2['isTerminal']==True,'WaterShedAreaKM2']), s=1)
# plt.show()
# plt.close()

# joined2['ton_excess_N'] = joined2.apply(lambda x: x['LocalExcess_TN_rv'] > 0, axis=1)
# plt.scatter(joined2['ton_MAyld_method_2'], joined2['TNYield'],s=1, c=joined2['ton_excess_N'])
# plt.show()
# plt.close()

# joined2['in_excess_N'] = joined2.apply(lambda x: x['excess_N_percent'] > 0, axis=1)
# plt.scatter(joined2['MAL_N'] / joined2['watershed_area_ha'], joined2['Load_N'] / joined2['watershed_area_ha'],s=1, c=joined2['in_excess_N'])
# plt.show()

# plt.scatter(random_thousand['Load_N'] / (100 * random_thousand['WaterShedAreaKM2']), random_thousand['TNYield'], s=1)
# plt.savefig('plots/N_load_scatter.png')
# plt.close()


