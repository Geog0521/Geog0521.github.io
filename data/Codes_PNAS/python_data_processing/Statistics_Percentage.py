from osgeo import gdal,gdalconst
import numpy as np
from tqdm import tqdm

def get_lc_30m(pixel_x_1km,pixel_y_1km,ds_1000m,ds_30m):
    """
    the function aims to get land cover classes at 30m
    :param pixel_x_1km: the pixel x coordinate at 1km
    :param pixel_y_1km: the pixel y coordinate at 1km
    :param ds_1000m: dataset of 1km land cover map opened by gdal
    :param ds_30m:   dataset of 30m land cover map opened by gdal
    :return:
    """

    geotransform_1000m = ds_1000m.GetGeoTransform()

    geotransform_30m = ds_30m.GetGeoTransform()

    #1km pixel coordinate
    pixel_x_1000m = pixel_x_1km #e.g.,300
    pixel_y_1000m = pixel_y_1km #e.g.,300

    # get geographical coordinate
    lon_1000m = geotransform_1000m[0] + pixel_x_1000m * geotransform_1000m[1]
    lat_1000m = geotransform_1000m[3] + pixel_y_1000m * geotransform_1000m[5]

    # get pixel coordinate
    pixel_x_30m = int((lon_1000m - geotransform_30m[0]) / geotransform_30m[1])
    pixel_y_30m = int((lat_1000m - geotransform_30m[3]) / geotransform_30m[5])

    resolution_ratio = geotransform_1000m[1] / geotransform_30m[1]
    resolution_ratio
    #the size of rectangle box
    rectangle_width = int(resolution_ratio)
    rectangle_height = int(resolution_ratio)

    #the top-left coordinate of rectangle box
    rectangle_x_start = int(pixel_x_30m - (rectangle_width // 2))
    rectangle_y_start = int(pixel_y_30m - (rectangle_height // 2))

    # print("left-top x:{}".format(rectangle_x_start))
    # print("left-top y:{}".format(rectangle_y_start))
    # confirm the geographic extent
    # 获取tif文件的宽度和高度
    tif_width = ds_30m.RasterXSize
    tif_height = ds_30m.RasterYSize

    # 确保x_start和y_start不超出tif的范围
    x_start = max(0, min(rectangle_x_start, tif_width - rectangle_width))
    y_start = max(0, min(rectangle_y_start, tif_height - rectangle_height))

    band = ds_30m.GetRasterBand(1) #get the first band
    try:
        landcover_values_30m = band.ReadAsArray(
            x_start, y_start, rectangle_width, rectangle_height)
    except Exception as e:
        print(e)
        exit()

    # if rectangle_width>=ds_30m.RasterXSize or rectangle_height>=ds_30m.RasterYSize:
    #     landcover_values_30m = None
    #
    # if rectangle_x_start >= ds_30m.RasterXSize or rectangle_y_start >= ds_30m.RasterYSize:
    #     landcover_values_30m = None
    #
    # else:
    #     # get all pixel values within the rectangle box
    #     try:
    #         print("get 30m pixels!")
    #
    #     except Exception as e:
    #         print(e)
    #         exit()

    return landcover_values_30m

def lc_cover_statistics(tif_30m, tif_1k, lc_list,year):
    """
    the function is developed to get the land cover percentage at 30m
    :param tif_30m: 30m land cover map in GeoTIFF file path
    :param tif_1k: this tif is the one upscaled from "tif_30m" with gdal.wrap
    :param lc_list: a list records the code for land cover class
    :param year_lis: store years, []
    :return:
    """
    #the 1km land cover map
    output_1k_mlc = tif_1k[:-4] + '_lc_perc.tif'
    import os
    if os.path.exists(output_1k_mlc):
        os.remove(output_1k_mlc)
    else:
        pass

    ds_30m = gdal.Open(tif_30m,gdalconst.GA_ReadOnly)

    # Open the input dataset
    ds_1k = gdal.Open(tif_1k, gdalconst.GA_ReadOnly)

    if ds_1k is None:
        print(f"Could not open input file: {tif_1k}")
        exit(1)

    xsize_1k = ds_1k.RasterXSize

    ysize_1k = ds_1k.RasterYSize

    geotransform = ds_1k.GetGeoTransform()

    # Create an output dataset with the desired pixel size
    driver = gdal.GetDriverByName('GTiff')

    # represent the number of land cover classes
    num_specific_classes = len(lc_list)

    #create a output dataset of specificed size and the number of bands
    output_ds = driver.Create(output_1k_mlc, xsize_1k, ysize_1k,  num_specific_classes, gdalconst.GDT_Float32)

    # Set the projection and geotransform for the output dataset
    output_ds.SetGeoTransform((geotransform[0], 1000, 0, geotransform[3], 0, -1000))

    output_ds.SetProjection(ds_1k.GetProjection())

    if output_ds is None:
        print(f"Could not create output file: {output_1k_mlc}")
        exit(1)

    # land cover map template at 1km resolution
    ds_1k = gdal.Open(tif_1k)

    # # land cover map at 30m resolution
    # ds_30m = gdal.Open(tif_30m)

    #statistics of land cover classes
    for lc_index in tqdm(range(len(lc_list)), desc="Outer Loop"):

        class_perc_arr = np.zeros((ysize_1k, xsize_1k), dtype=np.float32)

        for pixel_x_1k in range(ysize_1k): #row

            for pixel_y_1k in (range(xsize_1k)): #column
                                    #colunmn    #row
                window = get_lc_30m(pixel_y_1k, pixel_x_1k, ds_1k, ds_30m)
                window
                # Calculate the percentage for the specific land cover class
                percentage = 0
                if window is None:
                    percentage = 0
                    continue
                else:
                    count = np.sum(window == lc_list[lc_index])
                    total_pixels = len(window.flatten())
                    percentage = (count / total_pixels)
                #print("land cover percentage is :{}".format(percentage))
                try:
                    class_perc_arr[pixel_x_1k, pixel_y_1k] = percentage
                except Exception as e:
                    print(e)
                #print("percentage:{}".format(percentage))
        try:
            # import matplotlib.pyplot as plt
            # # plot the input image
            # plt.imshow(class_perc_arr)
            # plt.show()
            output_ds
            class_perc_arr = np.reshape(class_perc_arr,(ysize_1k,xsize_1k))

            output_ds.GetRasterBand(lc_index+1).WriteArray(class_perc_arr)

            #write band description
            output_ds.GetRasterBand(lc_index+1).SetDescription(str(year)+"_"+str(lc_list[lc_index]))
            print("\n percentage :{} information have been recorded! \n".format(lc_list[lc_index]))

        except Exception as e:
            print(e)
            exit()

    #delete memory
    ds_1k = None
    output_ds = None
    ds_30m = None

    return True

if __name__ == '__main__':
    # GeoTIFF file path
    # tif_30m = "output_lambert.tif"

    # this tif is the one upscaled from "tif_30m" with gdal.wrap
    # tif_1k = "1k_output_lambert.tif"
    #it is just a template, so we can use it repeatedly

    #a list records the code for land cover class
    # zenodo.org/record/8239305 for more details
    # Select NYC, Boston, San Diego and Dallas
    # 10 represents the Rainfed cropland
    # 120 denotes the Shurbland
    lc_cover_ls = [10, 11, 12, 20, 51, 52, 61, 62, 71, 72, 81, 82, 91, 92, 120, 121, 122, 130,
              140, 150, 152, 153, 181, 182, 183, 184, 185, 186, 187, 190, 200, 201, 202, 210, 0, 250]

    # years_lis = [i for i in range(21)]
    # GFLC30m_dire = "/Volumes/T7 Shield/gee_single/urban_GLFC30m"
    # City_list = ["Boston", "Dallas", "LA_SD"]

    GFLC30m_dire = "/Users/patrickchan/Downloads/gee_single/urban_GLFC30m"
    City_list = ["FL", "GA"]

    year_lis = [year for year in range(2000, 2020, 1)]  # 2000-2019
    for city in City_list:
        for year in year_lis:
            tif_30m = GFLC30m_dire +'/'+city+"/"+city+"_"+str(year)+"_lambert.tif"
            print(f"processing {tif_30m}")
            tif_lambert1k = GFLC30m_dire +'/'+city+"/"+city+"_"+str(year)+"_lambert_p1k.tif"
            print(f"processing {tif_lambert1k}")
            lc_cover_statistics(tif_30m, tif_lambert1k, lc_cover_ls,year)