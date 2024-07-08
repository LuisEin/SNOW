#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 12:52:47 2024

@author: le

This script takes Sentinel1 SAR Wet Snow data .zip files and
Unzips them
Opens the geotiff  and clips it to the AOI
Calculates the mean and the sum of wetsnow pixels
Saves the clipped file
Saves the mean and sum into a pandas_df with date_index

Depending on if the AOI is defined as a shapefile or ASCII file it needs to be 
changed in the bounds_aoi variable

"""
import os, glob, sys, zipfile, math
from osgeo import gdal, gdalconst, osr, ogr
from datetime import datetime as dt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


data_folder = "/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/Sentinel_Data/SWS/SWS_only_two_dates"
# 04_AOI_shapefile_Zugspitze_Watershed or 03_AOI_shp_zugspitze_reproj_for_code
aoi_path = '/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/Shapefiles/shapefile_Zugspitze/03_AOI_shp_zugspitze_reproj_for_code/AOI_zugspitze_reproj_32632.shp' #/shapefile_new_approach/mask_catchments_32632.asc
# This is for an ASCII file
mask_file_path = "/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/Shapefiles/shapefile_Zugspitze/04_AOI_shapefile_Zugspitze_Watershed/shapefile_new_approach/mask_catchments_32632.asc"
local_temp = "/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/Sentinel_Data/code/temp_new"
log_path = "/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/Sentinel_Data/code/logfile.txt"
output_folder = '/home/luis/Data/04_Uni/03_Master_Thesis/SNOW/02_data/Sentinel_Data/SWS/SWS_AOI_cropped/SWS_ASCII'



def debug_log(message):
    print(message)
    with open(log_path, "a") as log_file:
        log_file.write(f"{dt.now()}: {message}\n")

# New function to read the ASCII grid file
def read_ascii_grid(filename):
    with open(filename, 'r') as file:
        header = [file.readline() for _ in range(6)]
    ncols = int(header[0].split()[1])
    nrows = int(header[1].split()[1])
    xllcorner = float(header[2].split()[1])
    yllcorner = float(header[3].split()[1])
    cellsize = float(header[4].split()[1])
    nodata_value = float(header[5].split()[1])
    data = np.loadtxt(filename, skiprows=6)
    return ncols, nrows, xllcorner, yllcorner, cellsize, nodata_value, data

# New function to get bounds from the mask file
def get_bounds_from_mask(mask_file_path):
    ncols, nrows, xllcorner, yllcorner, cellsize, nodata_value, data = read_ascii_grid(mask_file_path)
    xmin = xllcorner
    xmax = xllcorner + ncols * cellsize
    ymin = yllcorner
    ymax = yllcorner + nrows * cellsize
    return [xmin, xmax, ymin, ymax]

def df_from_dir(directory):
    filelist = glob.glob("{}{}*.zip".format(directory, os.path.sep))

    if not all([os.path.basename(file).startswith("SWS") for file in filelist]):
        print("Not all Zipfiles seem to be FSC products! Check Data Folder! \nAborting Script execution...")
        sys.exit()

    filelist = [os.path.basename(file).split(".")[0] for file in glob.glob("{}{}*.zip".format(directory, os.path.sep))]
    sens_dates = [dt.strptime(file.split("_")[1], "%Y%m%dT%H%M%S") for file in filelist]
    sens_dates.sort()
    tiles = [file.split("_")[3] for file in filelist]

    scenes_df = pd.DataFrame({"tilecode": tiles, "sensdatetime": sens_dates, "filename": filelist})
    scenes_df["sensdate"] = scenes_df["sensdatetime"].round("h")

    return scenes_df

def writeLog(logfilepath, string, verbose=True):
    timestring = str(dt.now())[0:19]

    if os.path.isfile(logfilepath):
        with open(logfilepath, "a") as writefile:
            writefile.write("\n{}: {}".format(timestring, string))
        writefile.close()
    else:
        with open(logfilepath, "w") as writefile:
            writefile.write("\n{}: {}".format(timestring, string))
        writefile.close()
    if verbose:
        print(string)

def clipArray(src_data, params, bound, offsetx=0, offsety=0):
    res = params[2]
    srcMinX = params[0]
    srcMaxY = params[1]
    minX = bound[0]
    maxX = bound[1]
    minY = bound[2]
    maxY = bound[3]
    rowMin = int(math.ceil((srcMaxY - minY) / res))
    rowMax = int(round((srcMaxY - maxY) / res))
    colMin = int(round((minX - srcMinX) / res))
    colMax = int(math.ceil((maxX - srcMinX) / res))
    clip = src_data[rowMax + offsety:rowMin + offsety, colMin + offsetx:colMax + offsetx]
    minXclip = srcMinX + colMin * res
    maxYclip = srcMaxY - rowMax * res
    geoTrans = [minXclip, res, 0, maxYclip, 0, -res]
    
    debug_log(f"Clipping array with bounds: {bound}")
    debug_log(f"Clipped data shape: {clip.shape}")
    
    return clip, geoTrans

def getBounds_Shp(shapefile):
    shapeData = ogr.Open(shapefile, 0)
    layer = shapeData.GetLayer()
    bounds = layer.GetExtent()
    return bounds

def getBounds_Raster(ds):
    geoTrans = ds.GetGeoTransform()
    width = ds.RasterXSize
    height = ds.RasterYSize
    res = geoTrans[1]
    minX = geoTrans[0]
    maxX = minX + width * res
    maxY = geoTrans[3]
    minY = maxY - height * res
    bounds = [minX, maxX, minY, maxY]
    return bounds

def readRaster(filename):
    ds = gdal.Open(filename)
    band = ds.GetRasterBand(1)
    data = band.ReadAsArray()
    projection = osr.SpatialReference()
    projection.ImportFromWkt(ds.GetProjectionRef())
    geotrans = ds.GetGeoTransform()
    bounds = getBounds_Raster(ds)

    del ds, band
    return data, geotrans, projection, bounds

def getClipParams(ds):
    geoTrans_src = ds.GetGeoTransform()
    res = geoTrans_src[1]
    srcMinX = geoTrans_src[0]
    srcMaxY = geoTrans_src[3]
    params = [srcMinX, srcMaxY, res]
    return params

def write_grid(filename, data, geotrans, projection, driver="GTiff", dtype=None):
    debug_log(f"Writing grid to file: {filename}")
    debug_log(f"Data shape: {data.shape}, Geotransform: {geotrans}, Projection: {projection.ExportToWkt()}")

    if dtype is None:
        dtype = gdalconst.GDT_Float32

    driver = gdal.GetDriverByName(driver)
    if driver is None:
        debug_log("Driver GTiff is not available.")
        raise ValueError("Driver GTiff is not available.")
    
    debug_log(f"Creating file with dimensions: {data.shape[1]}, {data.shape[0]}")
    ds = driver.Create(filename, data.shape[1], data.shape[0], 1, dtype, ["COMPRESS=LZW"])
    if ds is None:
        debug_log("Failed to create the dataset.")
        raise ValueError("Failed to create the dataset.")
    
    ds.SetGeoTransform(geotrans)
    ds.SetProjection(projection.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(data)

    debug_log("Finished writing the grid.")
    del ds

def create_unique_filename(output_folder, base_name, ext):
    counter = 1
    unique_name = f"{base_name}.{ext}"
    while os.path.exists(os.path.join(output_folder, unique_name)):
        unique_name = f"{base_name}_{counter}.{ext}"
        counter += 1
    return os.path.join(output_folder, unique_name)

scenes_df = df_from_dir(data_folder)

if 'filename' not in scenes_df.columns:
    raise KeyError("The 'filename' column is missing from the DataFrame")

scenes_df['scene_date'] = scenes_df['filename'].apply(lambda x: dt.strptime(x.split("_")[1], "%Y%m%dT%H%M%S"))

if 'scene_date' not in scenes_df.columns:
    raise KeyError("The 'scene_date' column was not added to the DataFrame")

scenes_df = scenes_df.sort_values(by='scene_date')

days = scenes_df['scene_date'].dt.date.unique()

df_datestamp = pd.DataFrame(index=pd.Index(days, name='Date'))
df_datestamp["wetsnow_mean"] = np.nan
df_datestamp["wetsnow_sum"] = np.nan

bounds_aoi = getBounds_Shp(aoi_path)

for index, row in scenes_df.iterrows():
    scene_name = row["filename"]
    scene_date = row["scene_date"]

    debug_log(f"\nProcessing {scene_date} with scene {scene_name}\n")

    base_outfile_name = f"SWS_{scene_date.strftime('%Y_%m_%d_%H_%M_%S')}"
    outfile = create_unique_filename(output_folder, base_outfile_name, "tif")

    if os.path.isfile(outfile):
        debug_log(f"File {os.path.basename(outfile)} already exists, skipping processing!\n")
        continue

    zipfilepath = f"{data_folder}{os.path.sep}{scene_name}.zip"
    fileoutpath = f"{local_temp}{os.path.sep}{scene_name}.tif"

    try:
        with zipfile.ZipFile(zipfilepath) as z:
            with open(fileoutpath, "wb") as f:
                f.write(z.read(f"{scene_name}/{scene_name}_WSM.tif"))
        debug_log(f"Unzipped {zipfilepath} to {fileoutpath}\n")
    except Exception as e:
        writeLog(log_path, f"Couldnt unzip file {e} for following reason: {zipfilepath}")

    scenes = glob.glob(f"{local_temp}{os.path.sep}SWS*.tif")
    
    for scene in scenes:
        debug_log(f"Reading scene {scene}\n")
        rasterarray, scene_geotrans, projection, bounds = readRaster(scene)
        scene_params = [scene_geotrans[0], scene_geotrans[3], scene_geotrans[1]]
        data_aoi, geotrans_aoi = clipArray(rasterarray, scene_params, bounds_aoi)
        
        debug_log(f"Clip parameters: {scene_params}\nBounds: {bounds_aoi}\n")
        debug_log(f"Data AOI shape: {data_aoi.shape}\n")

        if np.all(data_aoi == 255):
            debug_log(f"Skipping and removing scene {scene} as it contains only NaN values.\n")
            continue

        data_aoi[data_aoi == 110] = 1
        data_aoi[data_aoi == 125] = 0
        data_aoi[data_aoi >= 200] = 255

        meanwetsnowarea = np.nanmean(data_aoi[data_aoi != 255])
        sumwetsnowpixels = np.nansum(data_aoi[data_aoi != 255])

        debug_log(f"Mean wetsnow part of scene is: {meanwetsnowarea}\n")
        debug_log(f"Wetsnow sum of scene is: {sumwetsnowpixels}\n")

        if not np.isnan(meanwetsnowarea) and not np.isnan(sumwetsnowpixels):
            debug_log(f"Writing grid to {outfile}\n")
            write_grid(outfile, data_aoi, geotrans_aoi, projection, dtype=gdal.GDT_Byte)
            df_datestamp.loc[scene_date.date(), 'wetsnow_mean'] = meanwetsnowarea
            df_datestamp.loc[scene_date.date(), 'wetsnow_sum'] = sumwetsnowpixels

    debug_log("Cleaning up temporary files\n")
    files = glob.glob(f"{local_temp}{os.path.sep}*")
    for f in files:
        os.remove(f)

plt.figure(figsize=(10, 6))
plt.plot(df_datestamp.index, df_datestamp['wetsnow_mean'], marker='o', linestyle='', color='blue')
plt.title('Mean proportion of Wetsnow Pixels Over Time')
plt.xlabel('Date')
plt.ylabel('Mean of Wetsnow Pixels')
plt.grid(True)
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(df_datestamp.index, df_datestamp['wetsnow_sum'], marker='o', linestyle='', color='blue')
plt.title('Total Sum of Wetsnow Pixels Over Time')
plt.xlabel('Date')
plt.ylabel('Sum of Wetsnow Pixels')
plt.grid(True)
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()
