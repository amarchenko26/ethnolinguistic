#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 14:51:31 2023

2) Show that for Thiessen polygons drawn around randomly placed points within 
continents (for the number of points in each continent, use the true number of groups), 
the empirical relationship between geographic variability at the country level and 
ethnic heterogeneity based on your Thiessen polygons does not hold. Produce the distribution of 
coefficient estimates based on 500 permutations of random points. 
Add the estimated coefficient based on the (1).

@author: anyamarchenko
"""

import os
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import random
from shapely.geometry import Point, Polygon
from scipy.spatial import Voronoi, voronoi_plot_2d
import subprocess # for sleep
import statsmodels.api as sm


def prevent_sleep():
    return subprocess.Popen(['caffeinate'])

def allow_sleep(process):
    process.terminate()

# Uncomment to prevent comp from sleeping
process = prevent_sleep()

# Set the base directory
base_directory = "/Users/anyamarchenko/Desktop/ethnolinguistic"
os.chdir(base_directory)

# Load the shapefiles
language_gdf = gpd.read_file('ethnologue/Ethnologue_16_shapefile/langa_no_overlap_biggest_clean.shp')

# =============================================================================
# # Re-project if the CRS is geographic
# if language_gdf.crs.is_geographic:
#     # Example: using a World Mercator projection
#     language_gdf = language_gdf.to_crs('EPSG:3395')
#     
# =============================================================================
    
languages_per_continent = language_gdf.groupby('CNT').size()


def generate_random_points(geometry, num_points):
    points = []
    minx, miny, maxx, maxy = geometry.bounds
    while len(points) < num_points:
        point = Point(random.uniform(minx, maxx), random.uniform(miny, maxy))
        if geometry.contains(point):
            points.append(point)
    return points


# Generate random points for each continent
random_points = {}
for continent, num_languages in languages_per_continent.items():
    # Get the combined geometry for all languages in this continent
    continent_geometry = language_gdf[language_gdf['CNT'] == continent].unary_union
    
    # Now generate random points within this geometry
    random_points[continent] = generate_random_points(continent_geometry, num_languages)


# Create a DataFrame to store results
polygon_languages_df = pd.DataFrame(columns=['Continent', 'Polygon', 'Num_Languages', 'Countries'])

# Generate Voronoi polygons and count languages
for continent, points in random_points.items():
    vor = Voronoi([point.coords[0] for point in points])
    polygons = [Polygon(vor.vertices[region]) for region in vor.regions if -1 not in region and region]

    for poly in polygons:
        contained_languages = language_gdf[language_gdf.intersects(poly)]
        countries = contained_languages['C1'].unique()  # Extract unique country codes
        new_row = pd.DataFrame([{
            'Continent': continent, 
            'Polygon': poly, 
            'Num_Languages': len(contained_languages),
            'Countries': ', '.join(countries)  # Join country codes as a string
        }])
        polygon_languages_df = pd.concat([polygon_languages_df, new_row], ignore_index=True)

# Display the DataFrame
print(polygon_languages_df)


### Clean up polygon_languages_df
# Define a function to randomly select a country from the string
def select_random_country(countries_str):
    countries_list = countries_str.split(', ')
    if countries_list:
        return random.choice(countries_list)
    return None

# Apply this function to the 'Countries' column and create a new column 'countryname'
polygon_languages_df['countryname'] = polygon_languages_df['Countries'].apply(select_random_country)

# Replace 'Russian Federation' with 'Russia' in the 'countryname' column
polygon_languages_df['countryname'] = polygon_languages_df['countryname'].replace('Russian Federation', 'Russia')
polygon_languages_df['countryname'] = polygon_languages_df['countryname'].replace('Viet Nam', 'Vietnam')
polygon_languages_df['countryname'] = polygon_languages_df['countryname'].replace('Iran', 'Iran, Islamic Rep.')
polygon_languages_df['countryname'] = polygon_languages_df['countryname'].replace('Egypt', 'Egypt, Arab Rep.')


### Merge with country agriculture data 
# Load the .dta file
dta_file_path = 'data/Tables1-3a.dta'  # Replace with the actual file path
country_data_df = pd.read_stata(dta_file_path)

# Select only the required columns from the country data
country_data_df = country_data_df[['countryname', 'sd_emeanclip', 'emeanclip', 'sdclimclip', 'sd_suitclip', 'abs_latclip']]

# Merge the data into polygon_languages_df
merged_df = polygon_languages_df.merge(country_data_df, on='countryname', how='left')

merged_df = merged_df.drop(columns=['Polygon', 'Countries'])
merged_df['Num_Languages'] = pd.to_numeric(merged_df['Num_Languages'], errors='coerce')


# Save the merged DataFrame to a new .dta file
merged_df.to_stata('data/thiessen_polygon.dta')


# Convert the relevant columns to numeric data types
merged_df['sd_emeanclip'] = pd.to_numeric(merged_df['sd_emeanclip'], errors='coerce')
merged_df['sdclimclip'] = pd.to_numeric(merged_df['sdclimclip'], errors='coerce')
merged_df['sd_suitclip'] = pd.to_numeric(merged_df['sd_suitclip'], errors='coerce')
merged_df['Num_Languages'] = pd.to_numeric(merged_df['Num_Languages'], errors='coerce')

# Display the merged DataFrame
print(merged_df)


merged_df = merged_df.dropna(subset=['Num_Languages', 'sd_emeanclip', 'sdclimclip', 'sd_suitclip'])

# Define the independent variables (predictors) and the dependent variable
X = merged_df[['sd_emeanclip', 'sdclimclip', 'sd_suitclip']]
y = merged_df['Num_Languages']

# Adding a constant to the model (for the intercept)
X = sm.add_constant(X)

# Create a model
model = sm.OLS(y, X)

# Fit the model
results = model.fit()

# Display the coefficient table
print(results.summary())




# =============================================================================
# # Create a DataFrame to store results
# polygon_languages_df = pd.DataFrame(columns=['Continent', 'Polygon', 'Num_Languages', 'Largest_Intersecting_Country'])
# 
# # Generate Voronoi polygons and count languages
# for continent, points in random_points.items():
#     vor = Voronoi([point.coords[0] for point in points])
#     polygons = [Polygon(vor.vertices[region]) for region in vor.regions if -1 not in region and region]
# 
#     # Create a GeoDataFrame for Voronoi polygons
#     voronoi_gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(polygons))
#     # Set the CRS for Voronoi polygons to match the language data
#     voronoi_gdf.crs = language_gdf.crs
# 
#     for poly in voronoi_gdf.geometry:
#         contained_languages = language_gdf[language_gdf.intersects(poly)].copy()
#         # Check if there are any intersecting languages
#         if contained_languages.empty:
#             print(f"No intersecting languages for polygon in {continent}")
#             continue
# 
#         # Calculate intersection area for each country
#         intersection_areas = contained_languages.intersection(poly).area
#         contained_languages['Intersection_Area'] = intersection_areas
#         # Find the country with the largest intersection
#         largest_country = contained_languages.loc[contained_languages['Intersection_Area'].idxmax(), 'C1']
# 
#         new_row = pd.DataFrame([{
#             'Continent': continent, 
#             'Polygon': poly, 
#             'Num_Languages': len(contained_languages),
#             'Largest_Intersecting_Country': largest_country
#         }])
#         polygon_languages_df = pd.concat([polygon_languages_df, new_row], ignore_index=True)
# 
# # Display the DataFrame
# print(polygon_languages_df)
# =============================================================================

# allow comp to sleep once code is done
allow_sleep(process)




