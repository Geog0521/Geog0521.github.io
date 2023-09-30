
#Combine four images in order to get image data for NYC
import os
from osgeo import gdal
import numpy as np
from tqdm import tqdm

def readTif(fileName):

    dataset = gdal.Open(fileName)
    if dataset == None:
        print("cannot open files {}".format(fileName))
        exit()
    else:
        return dataset

    return True

def writeTiff(im_data, im_geotrans, im_proj, path):
    im_bands = 0
    if 'int8' in im_data.dtype.name:
        datatype = gdal.GDT_Byte
    elif 'int16' in im_data.dtype.name:
        datatype = gdal.GDT_UInt16
    else:
        datatype = gdal.GDT_Float32
    if len(im_data.shape) == 3:
        im_bands, im_height, im_width = im_data.shape
    elif len(im_data.shape) == 2:
        im_data = np.array([im_data])
        im_bands, im_height, im_width = im_data.shape

    # create files
    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(path, int(im_width), int(im_height), int(im_bands), datatype)

    dataset.SetGeoTransform(im_geotrans)
    dataset.SetProjection(im_proj)
    for i in range(im_bands):
        dataset.GetRasterBand(i + 1).WriteArray(im_data[i])
    del dataset

    return True

def write_img_tif(tif_path,band_order_list,city):
    ##tif_path : image path
    ##band_order represent which band will be extracted
    dataset_img = readTif(tif_path)
    width = dataset_img.RasterXSize
    height = dataset_img.RasterYSize
    proj = dataset_img.GetProjection()
    geotrans = list(dataset_img.GetGeoTransform())

    file_path,file_name = os.path.split(tif_path)

    img_ds = dataset_img.ReadAsArray(0, 0, width, height)  #acquire data
    for band_order in band_order_list:
        saved_tif_dire = file_path + "/" + city
        saved_tif_dire
        if os.path.exists(saved_tif_dire):
            pass
        else:
            os.makedirs(saved_tif_dire)
            #/Users/patrickchan/Downloads/gee_single/urban_GLFC30m/Boston
            print(f"folder:'{saved_tif_dire}' has been created!")
        saved_tif = saved_tif_dire + "/" + city + "_" + str(2000+band_order) +".tif"
        # if os.path.exists(saved_tif):
        #     pass
        # else:
        #     os.remove(saved_tif)
        print(f"to save {saved_tif}")
        writeTiff(img_ds[band_order, :, :], geotrans, proj, saved_tif)  # band_order starts from 0

    dataset_img = None

def extract_band_year(GLFC30m_dire,GLFC30m_tif_lis,City_list):
    """
    This function aims to separate a multi-spectral image into individual bands(years)
    :param GLFC30m_file:
    :param ***:
    :return:
    """

    band_order_lis = [order for order in range(20)] #0-18
    for i in tqdm(range(len(GLFC30m_tif_lis)),desc="outer loop"):
        GLFC_tif_city = GLFC30m_dire+"/"+GLFC30m_tif_lis[i] #tif file
        print(f"It is processing'{GLFC_tif_city}'")
        write_img_tif(GLFC_tif_city, band_order_lis,City_list[i])

    print("finished!")

if __name__ == '__main__':
    # GFLC30m_dire = "/Users/patrickchan/Downloads/gee_single/urban_GLFC30m"
    # GLFC30m_tif_lis = ["GLC_FCS30D_20002022_W75N45_Annual.tif"
    #     ,"GLC_FCS30D_20002022_W100N35_Annual.tif",
    #             "GLC_FCS30D_20002022_W120N35_Annual.tif"]
    # City_list = ["Boston","Dallas","LA_SD"]
    # extract_band_year(GFLC30m_dire, GLFC30m_tif_lis,City_list)
    # print("processing data")

    #####
    GFLC30m_dire = "/Volumes/T7 Shield/back_up/urban_GLFC30m"
    GLFC30m_tif_lis = ["GLC_FCS30D_20002022_W85N30_Annual.tif"
        ,"GLC_FCS30D_20002022_W85N35_Annual.tif"]
    City_list = ["FL","GA"]
    extract_band_year(GFLC30m_dire, GLFC30m_tif_lis,City_list)
    print("processing data")

