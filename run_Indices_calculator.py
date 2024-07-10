#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 13:25:12 2024

This script loads 4Band PlanetScope Scenes
Extracts its bands 
Does Indices calculations with various bands
Exports the calculated index as .tif files
While also saving a copy of the clipped RGB image


@author: luis
"""

from osgeo import gdal
import numpy as np
import glob, os, shutil

def get_image_dimensions(file_path):
    dataset = gdal.Open(file_path)
    if not dataset:
        raise FileNotFoundError(f"Unable to open file: {file_path}")
    
    width = dataset.RasterXSize
    height = dataset.RasterYSize
    
    dataset = None
    return width, height

def mask_water(nir_band, green_band, threshold):
    """
    Masks water areas in the image based on a threshold calculated from green
    and nir bands. Based on Ji, L., Zhang, L., & Wylie, B. (2009). 
    Analysis of dynamic thresholds for the normalized difference water index. 
    Photogrammetric engineering & remote sensing, 75(11), 1307-1317.
    and Ji, L., Zhang, L., & Wylie, B. (2009). Analysis of dynamic thresholds for the normalized difference water index. Photogrammetric engineering & remote sensing, 75(11), 1307-1317.
    
    Parameters:
    red_band (numpy array): Array representing the red band of the image.
    threshold (float): Threshold value to identify water areas.
    
    Returns:
    numpy array: Masked array where water areas are set to NaN.
    
    Deactivated until further notice
    
    """
    water_mask = (green_band - nir_band)/(green_band + nir_band) > threshold
    return water_mask

def do_index_calculation(file, width, output_name, output_dir):
    dataset = gdal.Open(file)
    
    # Extract the date from the file name (assuming the date is part of the file name in the format YYYYMMDD)
    base_name = os.path.basename(file).split('.')[0]
    
    # Read the bands (assuming the bands are in the order: Blue, Green, Red, NIR)
    blue_band = dataset.GetRasterBand(1).ReadAsArray().astype(float)
    green_band = dataset.GetRasterBand(2).ReadAsArray().astype(float)
    red_band = dataset.GetRasterBand(3).ReadAsArray().astype(float)
    nir_band = dataset.GetRasterBand(4).ReadAsArray().astype(float)
    
    # # Apply water masking
    # water_mask = mask_water(nir_band, green_band, water_threshold)

    # Calculate indices based on output_name
    if output_name == "NDVI":
        index = (nir_band - red_band) / (nir_band + red_band)
    elif output_name == "BST":
        index = (nir_band - blue_band) / (nir_band + blue_band)
    elif output_name == "GST":
        index = (nir_band - green_band) / (nir_band + green_band)
    else:
        raise ValueError("Unknown index: {}".format(output_name))

    # Avoid division by zero
    index = np.where((nir_band + red_band) == 0, np.nan, index)
    
    # # Apply water mask
    # index[water_mask] = np.nan

    # Get georeference info
    geo_transform = dataset.GetGeoTransform()
    projection = dataset.GetProjection()

    # Create the output file for each tile
    driver = gdal.GetDriverByName('GTiff')
    out_file = os.path.join(output_dir, f'{base_name}_{output_name}_width_{width}px.tif')
    out_dataset = driver.Create(out_file, dataset.RasterXSize, dataset.RasterYSize, 1, gdal.GDT_Float32)

    # Set georeference info to the output file
    out_dataset.SetGeoTransform(geo_transform)
    out_dataset.SetProjection(projection)

    # Write the index band to the output file
    out_band = out_dataset.GetRasterBand(1)
    out_band.WriteArray(index)

    # Set NoData value
    out_band.SetNoDataValue(np.nan)

    # Flush data to disk
    out_band.FlushCache()
    out_dataset.FlushCache()

    # Clean up
    del dataset
    del out_dataset

    print(f"{output_name} calculation and export completed for {file}. Output saved as {out_file}")


# Define file paths and date pattern
# input_pattern = '/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/PlanetScope_Data/4_band/Data_March-June_23_Feb-March_24_psscene_analytic_sr_udm2/PSScene/*AnalyticMS_SR_clip.tif'  
input_pattern = "/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/PlanetScope_Data/merged_files/*_merged_masked.tif"
output_dir = '/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/PlanetScope_Data/code/temp/'
rgb_dir = os.path.join(output_dir, 'RGB')
ndvi_dir = os.path.join(output_dir, 'NDVI')
bst_dir = os.path.join(output_dir, 'BST')
gst_dir = os.path.join(output_dir, 'GST')


# List for the desired index calculations:
    # RGB
    # NDVI
    # BST - blue snow threshold
    # GST - green snow threshold
outputs = ["NDVI", "BST", "GST"]

# List for the desired output folders
out_dirs = [ndvi_dir, bst_dir, gst_dir]


# # Water masking threshold
# water_threshold = 0.0  # As 


# Create the output directories if they don't exist
os.makedirs(rgb_dir, exist_ok=True)
for dir in out_dirs:
    os.makedirs(dir, exist_ok=True)

# Main loop to iterate over the indices and the output directories

# filter the files, for only files covering the whole area of interest
    # Find all tiles for the specified pattern
    tile_files = glob.glob(input_pattern)
    
    # Filter files that cover the whole area and are wider than 200 pixels
    filtered_files = []
    for file in tile_files:
        width, height = get_image_dimensions(file)
        if width > 200 and height > 200:
            filtered_files.append((file, width))
            # Copy file to the RGB directory
            shutil.copy(file, rgb_dir)
    
    # Ensure we have at least one file to process
    if not filtered_files:
        raise ValueError("No suitable files found for the specified criteria.")
        
# Loop through each index and its corresponding directory
for output_name, output_dir in zip(outputs, out_dirs):
    for file, width in filtered_files:
        do_index_calculation(file, width, output_name, output_dir)