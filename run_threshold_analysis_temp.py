# -*- coding: utf-8 -*-
'''
Created on Tuesday Jul 09 09:13:12 2024

This file takes the BlueSnowThreshold calculated PlanetScope data and 
smooths it using a 1-d gaussian filter.
Then does Histograms and threshold calculations.

@author: luis
'''


from osgeo import gdal
from scipy.ndimage import gaussian_filter1d
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import os

def apply_gaussian_filter_and_generate_histogram(tiff_file, output_histogram_dir):
    # Open the TIFF file
    dataset = gdal.Open(tiff_file)
    if dataset is None:
        print(f"Unable to open {tiff_file}")
        return

    # Read the raster data
    band = dataset.GetRasterBand(1)
    array = band.ReadAsArray()

    # Apply Gaussian filter
    gau_array = gaussian_filter1d(array, 3)

    # Extract the date from the file name
    base_name = os.path.basename(tiff_file)
    date_str = base_name.split("_")[0]
    date = datetime.strptime(date_str, "%Y%m%d").strftime("%Y-%m-%d")

    # Plot Histogram
    plt.hist(gau_array.flatten(), bins=50)
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of BST from {date}')
    
    # Create the output directory if it doesn't exist
    date_dir = os.path.join(output_histogram_dir, date)
    if not os.path.exists(date_dir):
        os.makedirs(date_dir)
    
    # Construct the output file path
    output_file = os.path.join(date_dir, f'{date}_BST_hist.png')
    
    # Save the histogram
    plt.savefig(output_file)
    plt.close()
    print(f'Histogram saved to {output_file}')

    # Close the dataset
    dataset = None

def process_histograms(input_dir, output_histogram_dir):
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith('.tif'):
                tiff_file = os.path.join(root, file)
                apply_gaussian_filter_and_generate_histogram(tiff_file, output_histogram_dir)

if __name__ == "__main__":
    # Define the input directory containing the cropped TIFF files
    input_dir = "/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/PlanetScope_Data/Indices/merged_tiles/BST"

    # Define the output directory where histograms will be saved
    output_histogram_dir = "/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/PlanetScope_Data/Indices/BST_Histograms"

    # Process the histograms
    process_histograms(input_dir, output_histogram_dir)
