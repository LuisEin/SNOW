# -*- coding: utf-8 -*-
'''
Created on Tuesday Jul 09 17:13:12 2024

This file takes PlanetScope Scenes by the same date and puts the single tiles 
together.
Saves them as onetif file containing the whole AOI.

@author: luis
'''

import os
import glob
from osgeo import gdal, ogr
import numpy as np

# Define paths
input_pattern = '/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/PlanetScope_Data/4_band/Data_March-June_23_Feb-March_24_psscene_analytic_sr_udm2/PSScene/*AnalyticMS_SR_clip.tif'
output_dir = '/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/PlanetScope_Data/merged_files'
shapefile_path = '/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/Shapefiles/shapefile_Zugspitze/03_AOI_shp_zugspitze_reproj_for_code/AOI_zugspitze_reproj_32632.shp'  

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Collect files matching the pattern
files = glob.glob(input_pattern)

# Group files by date
files_by_date = {}
for file in files:
    basename = os.path.basename(file)
    date_str = basename[:8]
    if date_str not in files_by_date:
        files_by_date[date_str] = []
    files_by_date[date_str].append(file)

# Function to read shapefile
def read_shapefile(shapefile_path):
    shapefile = ogr.Open(shapefile_path)
    layer = shapefile.GetLayer()
    geom = layer.GetFeature(0).GetGeometryRef()
    return geom

# Read AOI from shapefile
aoi_geom = read_shapefile(shapefile_path)

# Process each date group
for date_str, file_list in files_by_date.items():
    print(f"Processing date: {date_str}")

    # List of datasets
    datasets = [gdal.Open(file) for file in file_list]

    # Get geotransform and projection from first dataset
    geotransform = datasets[0].GetGeoTransform()
    projection = datasets[0].GetProjection()

    # Merge datasets
    vrt_options = gdal.BuildVRTOptions()
    vrt = gdal.BuildVRT('/vsimem/temp.vrt', datasets, options=vrt_options)
    mosaic = gdal.Translate('', vrt, format='MEM')

    # Write merged file
    output_path = os.path.join(output_dir, f"{date_str}_merged.tif")
    gdal.Translate(output_path, mosaic)

    # Mask the merged file with the shapefile AOI
    output_masked_path = os.path.join(output_dir, f"{date_str}_merged_masked.tif")
    gdal.Warp(output_masked_path, output_path, cutlineDSName=shapefile_path, cropToCutline=True, dstNodata=np.nan)

    # Close datasets
    for dataset in datasets:
        dataset = None

print("Processing complete.")
