from osgeo import gdal, osr, gdalconst
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
import pyproj
from osgeo import gdal
import numpy as np
import os

def check_tif(tif_file):

    if os.path.exists(tif_file):
        os.remove(tif_file)
    else:
        pass

def set_projetion_pixel_size(input_tif_file):

    # Input and output file paths
    input_path = input_tif_file

    output_path = input_tif_file[:-4]+"_"+'lambert.tif'

    print(f" it is processing {output_path}")
    output_path
    try:
        check_tif(output_path)
    except Exception as e:
        print(e)

    # Define the source and target projections
    # source_crs = rasterio.crs.CRS.from_epsg(4326)  # Assuming input is in WGS84

    # Define the Proj4 string for EPSG:102009
    # https://epsg.io/102009
    target_proj4 = (
        '+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'
    )

    # #NA_Lambert_Conformal_Conic_2SP
    # target_proj4 = (
    #     '+proj=lcc +lat_1=60 +lat_2=25 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'
    # )

    # Create a CRS object from the Proj4 string
    target_crs = pyproj.CRS.from_string(target_proj4)

    # Open the input file
    with rasterio.open(input_path) as src:
        # Calculate the transform and dimensions for the target CRS
        transform, width, height = calculate_default_transform(src.crs, target_crs, src.width, src.height, *src.bounds)

        # Create the output file with the target CRS and dimensions
        with rasterio.open(output_path, 'w', driver='GTiff', crs=target_crs, transform=transform, width=width,
                           height=height,count=1,dtype=src.dtypes[0]) as dst:

            # Reproject the data
            reproject(
                source = rasterio.band(src, 1),  # Assuming you're working with the first band
                destination = rasterio.band(dst, 1),  # Assuming you're writing to the first band
                src_transform = src.transform,
                src_crs = src.crs,
                dst_transform = transform,
                dst_crs = target_crs,
                resampling = Resampling.nearest  # You can choose a different resampling method
            )
            print("reprojection is done!")

    return output_path

def resampling_to_1k(projected_tif):
    """
    the function that resampling a tif to 1k resolution
    :param in_tif: the input
    :param out_tif: the output
    :return:
    """
    # in_tif = "output_lambert.tif"
    out_tif = projected_tif [:-4]+"_p1k.tif"

    try:
        check_tif(out_tif)
    except Exception as e:
        print(e)

    x_res = 1000

    y_res = 1000

    resample_alg = gdal.GRA_NearestNeighbour

    resamp_resu = gdal.Warp(out_tif, projected_tif, xRes=x_res, yRes=y_res, resampleAlg=resample_alg)

    print("resampling to 1k resolution sucessfully")

    return True

def repro_resam_tif():
    GFLC30m_dire = "/Users/patrickchan/Downloads/gee_single/urban_GLFC30m"
    City_list = ["Boston", "Dallas", "LA_SD"]

    GFLC30m_dire = "/Users/patrickchan/Downloads/gee_single/urban_GLFC30m"
    City_list = ["FL", "GA"]

    year_lis = [year for year in range(2000, 2020, 1)]  # 2000-2019
    for city in City_list:
        for year in year_lis:
            input_tif = GFLC30m_dire +'/'+city+"/"+city+"_"+str(year)+".tif"
            print("Projection")
            projected_tif = set_projetion_pixel_size(input_tif)
            print("resampling")
            resampling_to_1k(projected_tif)
            print("finished!")
    return True
if __name__ == '__main__':

    # GeoTIFF file path
    repro_resam_tif()
