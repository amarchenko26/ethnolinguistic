#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 17:38:17 2023

@author: anyamarchenko
"""

import os
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.crs as ccrs


# Set the base directory
base_directory = "/Users/anyamarchenko/Desktop/ethnolinguistic"
os.chdir(base_directory)

# Load the shapefiles
virtual_country_gdf = gpd.read_file('data/Virtual_country/virtual_cntrygrid.shp')
language_gdf = gpd.read_file('ethnologue/Ethnologue_16_shapefile/langa_no_overlap_biggest_clean.shp')

# Drop langauges with fewer 1000 speakers
language_gdf = language_gdf[language_gdf['LMP_POP1'] >= 1000]

# Drop virtual countries with allcrops == 0 (aka no land quality data)
virtual_country_gdf = virtual_country_gdf[virtual_country_gdf['allcrops'] != 0]

# Plotting
fig, ax = plt.subplots(figsize=(10, 10))

# =============================================================================
# # Plot virtual countries
# virtual_country_gdf.plot(ax=ax, color='white', edgecolor='black')
# 
# # Plot languages
# language_gdf.plot(ax=ax, color='red', alpha=0.5)
# 
# plt.show()
# 
# =============================================================================

# Spatial join to find the intersection of languages and virtual countries
intersection_gdf = gpd.sjoin(virtual_country_gdf, language_gdf, how='inner', predicate='intersects')

# Calculate the number of languages per virtual country
languages_per_country = intersection_gdf.groupby('uniq_cnt25').size().reset_index(name='num_languages')

# Export the data to a Stata .dta file
languages_per_country.to_stata('data/languages_per_virtual_country.dta')

# Load the Stata file
table4_7b_df = pd.read_stata('data/Tables4-7b.dta')

# Merge the DataFrames on 'uniq_cnt25'
merged_df = pd.merge(table4_7b_df, languages_per_country, on='uniq_cnt25', how='left')

# Save the merged DataFrame to a new .dta file
merged_df.to_stata('data/Tables4-7b_merged.dta')