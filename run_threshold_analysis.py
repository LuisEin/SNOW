# -*- coding: utf-8 -*-
'''
Created on Tuesday Jul 09 09:13:12 2024

This file takes the BlueSnowThreshold calculated PlanetScope data and 
smooths it using a 1-d gaussian filter.
Then does Histograms and threshold calculations.

@author: luis
'''

from osgeo import gdal, ogr
from scipy.ndimage import gaussian_filter1d
from scipy.stats import gaussian_kde
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import glob, os, shutil




def get_aoi_bounds(shapefile):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    datasource = driver.Open(shapefile, 0)  # 0 means read-only. 1 means writeable.
    if not datasource:
        print(f"Unable to open {shapefile}")
        return None

    layer = datasource.GetLayer()
    envelope = layer.GetExtent()  # Returns (xmin, xmax, ymin, ymax)

    
    return (envelope[0], envelope[2], envelope[1], envelope[3])  # (xmin, ymin, xmax, ymax)


def get_bounds(shapefile):
    shapeData = ogr.Open(shapefile, 0)
    layer = shapeData.GetLayer()
    bounds = layer.GetExtent()
    return bounds



def check_and_crop_tiff(tiff_file, aoi_bounds, output_dir):
    # Open the TIFF file
    dataset = gdal.Open(tiff_file)
    if dataset is None:
        print(f"Unable to open {tiff_file}")
        return

    # Get georeference info
    transform = dataset.GetGeoTransform()
    pixel_width = transform[1]
    pixel_height = transform[5]
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize

    # Get the bounding box of the raster
    raster_bounds = (
        transform[0],  # xmin
        transform[3] + (rows * pixel_height),  # ymin
        transform[0] + (cols * pixel_width),  # xmax
        transform[3]  # ymax
    )

    # Check if AOI is covered
    if (aoi_bounds[0] < raster_bounds[2] and aoi_bounds[2] > raster_bounds[0] and
        aoi_bounds[1] < raster_bounds[3] and aoi_bounds[3] > raster_bounds[1]):
        
        # Calculate the pixel coordinates of the AOI
        x_offset = int((aoi_bounds[0] - transform[0]) / pixel_width)
        y_offset = int((aoi_bounds[3] - transform[3]) / pixel_height)
        x_size = int((aoi_bounds[2] - aoi_bounds[0]) / pixel_width)
        y_size = int((aoi_bounds[3] - aoi_bounds[1]) / abs(pixel_height))

        # Ensure offsets and sizes are within the raster bounds
        x_offset = max(x_offset, 0)
        y_offset = max(y_offset, 0)
        x_size = min(x_size, cols - x_offset)
        y_size = min(y_size, rows - y_offset)

        # Read the data from the specified window
        band = dataset.GetRasterBand(1)
        data = band.ReadAsArray(x_offset, y_offset, x_size, y_size)

        # Define the new geotransform for the cropped raster
        new_transform = (
            transform[0] + x_offset * pixel_width,
            pixel_width,
            0.0,
            transform[3] + y_offset * pixel_height,
            0.0,
            pixel_height
        )

        # Create the output directory if it doesn't exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Create the output file
        output_file = os.path.join(output_dir, os.path.basename(tiff_file))
        driver = gdal.GetDriverByName('GTiff')
        out_dataset = driver.Create(output_file, x_size, y_size, 1, band.DataType)
        out_dataset.SetGeoTransform(new_transform)
        out_dataset.SetProjection(dataset.GetProjection())
        out_band = out_dataset.GetRasterBand(1)
        out_band.WriteArray(data)
        out_band.FlushCache()

        print(f"Cropped file saved to {output_file}")
    else:
        print(f"AOI not covered by {tiff_file}")

    # Close the dataset
    dataset = None

def process_tiff_files(input_dir, aoi_bounds, output_dir):
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith('.tif'):
                tiff_file = os.path.join(root, file)
                check_and_crop_tiff(tiff_file, aoi_bounds, output_dir)

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


## Define paths and run functions

if __name__ == "__main__": # Wenn dieses Modul in ein anderes Skript importiert wird, wird __name__ umbenannt z.B. in den skript Namen
    # Path to the AOI shapefile
    shapefile = "/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/Shapefiles/shapefile_Zugspitze/03_AOI_shp_zugspitze_reproj_for_code/AOI_zugspitze_reproj_32632.shp"

    # Get the AOI bounds from the shapefile
    aoi_bounds = get_bounds(shapefile)
    if aoi_bounds is None:
        print("Failed to get AOI bounds")
        exit(1)
        
    # Define the input directory containing the TIFF files
    input_dir = "/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/PlanetScope_Data/Indices/merged_tiles/BST"

    # Define the output directory where cropped TIFF files will be saved
    output_dir = "/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/PlanetScope_Data/Indices/BST/AOI_cropped_files"

    # Process the TIFF files
    process_tiff_files(input_dir, aoi_bounds, output_dir)

    # Define the input directory containing the cropped TIFF files
    input_dir = output_dir

    # Define the output directory where histograms will be saved
    output_histogram_dir = "/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/PlanetScope_Data/Indices/BST_Histograms"

    # Process the histograms
    process_histograms(input_dir, output_histogram_dir)





# Define file paths and date pattern
input_pattern = '/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/PlanetScope_Data/Indices/BST/20240229_101353_47_2475_3B_AnalyticMS_SR_clip_BST_width_6679px.tif' 
# /home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/PlanetScope_Data/Indices/merged_tiles/BST/20240229_merged_masked_BST_width_3332px.tif
# Adjust the pattern to match your tiles
output_dir = '/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/PlanetScope_Data/code/temp/'

## Load image
# Using gdal
dataset = gdal.Open(input_pattern, gdal.GA_ReadOnly)
for x in range(1, dataset.RasterCount + 1):
    band = dataset.GetRasterBand(x)
    array = band.ReadAsArray()

## Filter 
gau_array = gaussian_filter1d(array, 3)

# check for no data values
nodata_value = band.GetNoDataValue()
print(f"Nodata-Wert: {nodata_value}")



## Plot Histogram 
plt.hist(gau_array)
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram of BST_merged from the 29 Feb 2024') #{input_pattern[75:83]}
plt.show()

### Plot probability density function
# Erst abflachen
flat_data = gau_array.ravel()
density = gaussian_kde(flat_data)
x = np.linspace(min(flat_data), max(flat_data), 1000)
y = density(x)
# Probability Density Plot erstellen
plt.plot(x, y, label='Density')
plt.fill_between(x, y, alpha=0.5)
plt.xlabel('BST')
plt.ylabel('Density')
plt.title('Probability Density Plot 18 Feb 2024')
plt.legend()
plt.show()

### Schnellerer approach?

# Daten binning
counts, bin_edges = np.histogram(gau_array, bins=100)

# Bin-Mitten berechnen
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# Bar-Plot erstellen
plt.bar(bin_centers, counts, width=bin_edges[1] - bin_edges[0])
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Binned Data Bar Plot BST 29 Feb 2024')
plt.show()
