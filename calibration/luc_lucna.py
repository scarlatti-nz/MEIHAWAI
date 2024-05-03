import geopandas
from numpy import pi
import pyogrio
import pandas as pd

# cross tabulation of land use, land use capability, and slope based on the LUCAS dataset, used for sense checking.


# data files not included in the repo, available for download from:
# lucas - https://data.mfe.govt.nz/layer/52375-lucas-nz-land-use-map-1990-2008-2012-2016-v011/
# lri - https://lris.scinfo.org.nz/layer/48076-nzlri-land-use-capability-2021/

gsd1 = pyogrio.read_dataframe("data/lucas-nz-land-use-map-1990-2008-2012-2016-v011.shp")
gsd2 = pyogrio.read_dataframe("data/nzlri-land-use-capability-2021.shp")

merged = geopandas.overlay(gsd1, gsd2, how='intersection')
merged['area'] = merged['geometry'].apply(lambda x: x.area)
merged = merged.drop(columns=['geometry'])

slopedf = merged.groupby(['slope', 'LUCNA_2016', 'LUCID_2016', 'SUBID_2016']).sum().reset_index()
lucdf = merged.groupby(['luc', 'LUCNA_2016', 'LUCID_2016', 'SUBID_2016']).sum().reset_index()
bothdf = merged.groupby(['slope', 'luc', 'LUCNA_2016', 'LUCID_2016', 'SUBID_2016']).sum().reset_index()


slopedf = slopedf[slopedf['LUCID_2016'].isin([72,73,75,76,77,78])]
slopedf['slope'] = slopedf['slope'].apply(lambda x: x[0])
pivoted = slopedf.pivot_table(index='slope', columns=['LUCNA_2016','SUBID_2016'], values='area', aggfunc='sum')
pivoted.to_csv('calibration/slope_lucna.csv')

lucdf = lucdf[lucdf['LUCID_2016'].isin([72,73,75,76,77,78])]
lucdf['luc'] = lucdf['luc'].apply(lambda x: x[0])
pivot_luc = lucdf.pivot_table(index='luc', columns=['LUCNA_2016','SUBID_2016'], values='area', aggfunc='sum')
pivot_luc.to_csv('calibration/luc_lucna.csv')

bothdf = bothdf[bothdf['LUCID_2016'].isin([72,73,75,76,77,78])]
bothdf['slope'] = bothdf['slope'].apply(lambda x: x[0])
bothdf['luc'] = bothdf['luc'].apply(lambda x: x[0])
pivot_both = bothdf.pivot_table(index=['slope', 'luc'], columns=['LUCNA_2016','SUBID_2016'], values='area', aggfunc='sum')
pivot_both.to_csv('calibration/both_lucna.csv')


# area_by_slope = slopedf[['area', 'slope']]
# slopedf_sum = slopedf.groupby('slope')['area'].sum()

# slope_classes = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
# slope_classified = {}

# for slope_class in slope_classes:
#     slope_classified[slope_class] = slopedf_sum[slopedf_sum.index.str.startswith(slope_class)].sum()


# area_by_LUC = lucdf[['LUCNA_2016', 'area']]
# lucna_classes = ['72 - Planted Forest - Pre 1990', '75 - Grassland - High producing', '76 - Grassland - Low producing', '78 - Cropland - Annual']

# filtered_area_by_LUC = area_by_LUC[area_by_LUC['LUCNA_2016'].isin(lucna_classes)]
# summed_areas = filtered_area_by_LUC.groupby('LUCNA_2016')['area'].sum().reset_index()

# summed_areas_df = geopandas.GeoDataFrame(summed_areas, columns=['LUCNA_2016', 'summed_area'])
# summed_slope_df = geopandas.GeoDataFrame(slope_classified.items(), columns=['slope', 'summed_area'])

# print(summed_areas_df)

# print(summed_slope_df)

# intersection_slope_LUNCA = merged.groupby(['LUCNA_2016', 'slope']).sum().reset_index()
# intersections = intersection_slope_LUNCA[['LUCNA_2016', 'slope', 'area']]






