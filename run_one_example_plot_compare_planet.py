# -*- coding: utf-8 -*-
'''
Created on Tuesday Jul 10 15:32:12 2024

This script takes an example of merged and clipped tiles from the run_merge_tiles.py
And finds a threshold for shady and non-shady areas

@author: luis
'''

from osgeo import gdal, ogr
from scipy.ndimage import gaussian_filter1d
from scipy.stats import gaussian_kde
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import glob, os, shutil



# Define file paths and date pattern
input_pattern = '/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/PlanetScope_Data/Indices/merged_tiles/BST/20240229_merged_masked_BST_width_3332px.tif' 
# RGB
input_pattern2 = '/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/PlanetScope_Data/Indices/merged_tiles/RGB/20240229_merged_masked.tif'
# /home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/PlanetScope_Data/Indices/merged_tiles/BST/20240229_merged_masked_BST_width_3332px.tif
# Adjust the pattern to match your tiles
output_dir = '/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/PlanetScope_Data/code/temp/'

## Load indexed image 
dataset = gdal.Open(input_pattern, gdal.GA_ReadOnly)
for x in range(1, dataset.RasterCount + 1):
    band = dataset.GetRasterBand(x)
    array = band.ReadAsArray()

## Filter 
gau_array = gaussian_filter1d(array, 3)

# Get unique values and their counts
unique_values, counts = np.unique(array, return_counts=True)

# Display the unique values and their counts
print("Unique values:", unique_values)
print("Counts:", counts)

# Plot the array using Matplotlib
plt.imshow(gau_array, cmap='gray')
plt.colorbar()
plt.title('BST gaussian filtered array 2024_03_29')
plt.show()

# Display Histogram
# Plot Histogram
plt.hist(array.flatten(), bins=50)
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram of BST array 2024_03_29')


# check for no data values
nodata_value = band.GetNoDataValue()
print(f"Nodata-Wert: {nodata_value}")


######## Load RGB image ####################################################### 
dataset = gdal.Open(input_pattern2, gdal.GA_ReadOnly)

# Read the bands
bands = []
for i in range(1, 5):
    band = dataset.GetRasterBand(i).ReadAsArray()
    bands.append(band)

# Close the dataset
dataset = None


######### Plot Histograms of all four channels ################################
# Plot histograms
band_names = ['Red', 'Green', 'Blue', 'NIR']
plt.figure(figsize=(12, 8))

# Adjust text size
plt.rcParams.update({'font.size': 10})

for i, band in enumerate(bands):
    plt.subplot(2, 2, i+1)
    plt.hist(band.flatten(), bins=256, color='gray', alpha=0.7)
    plt.suptitle('RGB PlanetScope Image 2024-03-29')
    plt.title(f'{band_names[i]} Band')
    plt.xlabel('Pixel Value')
    plt.ylabel('Frequency')

plt.tight_layout()
plt.show()


######### Plot a comparison of BSI with Coastal blue vs blue ##################