#this script aims to get statistics of ...
#table_list
#Tampa.csv, LA.csv , ,
#city, year, lc_id, percentage
import os
from osgeo import gdal
import numpy as np
import pandas as pd
from tqdm import tqdm

def readTif(fileName):

    dataset = gdal.Open(fileName)
    if dataset == None:
        print("cannot open files {}".format(fileName))
        exit()
    else:
        return dataset

    return True

def get_statisics(tif_path):

    """
    #get the statistics of each land cover percentage in one year
    :param tif_path: e.g., Boston_2000_lambert_p1k_lc_perc_2_final_2.tif
    :return: a list the records total sum of land cover percentage
    """

    #open a tif at one year
    dataset_img = readTif(tif_path)
    width = dataset_img.RasterXSize
    height = dataset_img.RasterYSize

    # file_path,file_name = os.path.split(tif_path)
    img_ds = dataset_img.ReadAsArray(0, 0, width, height)  # acquire data
    img_ds
    sum_band_sta_lis = []

    for band_index in range(36): #band_index indicate land cover
        # band_index = 0  # start from 1
        band_statistics = img_ds[band_index, :, :]
        # print(band_statistics.min())
        # print(band_statistics.max())
        band_statistics

        # 使用np.isnan函数标识NaN值
        nan_mask = np.isnan(band_statistics)

        non_nan_values = band_statistics[~nan_mask]
        non_nan_values
        sum_band_statistics = np.sum(non_nan_values) #get the total sum in a year
        sum_band_statistics
        #np.nansum()
        # if sum_band_statistics>=0:
            #print("percentage:{}".format(sum_band_statistics))
        sum_band_sta_lis.append(sum_band_statistics)

    dataset_img = None

    return sum_band_sta_lis

def save_stastic_city(saved_csv,lc_id_lis,lc_descri_lis,sum_year_statistics):

    #write statistics into a csv

    products_list = [lc_id_lis,lc_descri_lis,sum_year_statistics]
    products_list
    df = pd.DataFrame(products_list).transpose()

    df.columns=['lc_id', 'lc_descri','percentage']

    if os.path.exists(saved_csv):
        os.remove(saved_csv)
    else:
        pass
    try:
        df.to_csv(saved_csv, index=False)
        print(f"{saved_csv} has been created!")
    except Exception as e:
        print(f"{saved_csv} cannot be created")
        print(e)

    dataset_img = None

    return True

if __name__ == '__main__':

    lc_id_lis = [10, 11, 12, 20, 51, 52, 61, 62, 71, 72, 81, 82, 91, 92, 120, 121, 122, 130,
              140, 150, 152, 153, 181, 182, 183, 184, 185, 186, 187, 190, 200, 201, 202, 210, 0, 250]
    lc_descri_lis = ["rainfed cropland",
                     "herbaceous cover cropland",
                     "tree or shrub cover(orchard) cropland",
                     "irrigated cropland",
                     "open evergreen broadleaved forest",
                     "closed evergreen broadleaved forest",
                     "open deciduous broadleaved forest(0.15<fc<0.4)",
                     "open deciduous broadleaved forest(fc>0.4)",
                     "open evergreen needle-leaved forest(0.15<fc<0.4)",
                     "closed evergreen needle-leaved forest(fc>0.4)",
                     "open deciduous needle-leaved forest(0.15<fc<0.4)",
                     "closed needle-leaved forest(fc>0.4)",
                     "open mixed leaf forest(broadleaved and needle-leaved)",
                     "closed mixed leaf forest (broadleaved and needle-leaved)",
                     "shrubland",
                     "evergreen shrubland",
                     "deciduous shrubland",
                     "Grassland",
                     "lichens and mosses",
                     "sparse vegetation(fc<0.15)",
                     "sparse shrubland (fc<0.15)",
                     "sparse herbaceous (fc<0.15)",
                     "swamp",
                     "marsh",
                     "looded flat",
                     "saline",
                     "mangrove",
                     "salt marsh",
                     "tidal flat",
                     "impervious surfaces",
                     "bare areas",
                     "consolidated bare areas",
                     "unconsolidated bare areas",
                     "water body",
                     "filled value 0",
                     "filled value 250"] #"permanent ice and snow"
    #pay attention to FL
    print("the total number of land cover clases is {}".format(len(lc_descri_lis)))
    print("the total number of land cover clases is {}".format(len(lc_id_lis)))
    main_dire = "/Volumes/T7 Shield/back_up/urban_GLFC30m"
    city_lis = ["Boston","Dallas","LA_SD","GA", "FL", "FL"]
    real_city_lis = ["Boston","Dallas","LA","Augs", "Orlando", "Tampa"]
    for i in tqdm(range(6)):
        print("Processing city:{}".format(real_city_lis[i]))
        for year in range(2000,2020,1):
            if i==5:
                masked_GLFC_1k = main_dire + "/" + "FL2" + "/" + city_lis[i] + "_" + str(
                    year) + "_lambert_p1k_lc_perc_2_final_2.tif"
            else:
                masked_GLFC_1k = main_dire + "/" + city_lis[i] + "/" + city_lis[i] + "_" + str(
                    year) + "_lambert_p1k_lc_perc_2_final_2.tif"

            sum_year_statistics = get_statisics(masked_GLFC_1k)
            #generate .csv
            file_pa, file_name = os.path.split(masked_GLFC_1k)
            csv_saved_folder = "final_statistics"
            saved_csv_path = file_pa + "/" + csv_saved_folder
            if os.path.exists(saved_csv_path):
                pass
            else:
                os.makedirs(saved_csv_path)
            saved_csv = saved_csv_path +"/" + real_city_lis[i] + "_" + str(year) + "_lc_st.csv"
            save_stastic_city(saved_csv, lc_id_lis, lc_descri_lis, sum_year_statistics)



