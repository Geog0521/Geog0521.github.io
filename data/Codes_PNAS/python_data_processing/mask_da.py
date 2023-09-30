import xarray as xr
import rasterio as rio
import os
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal,gdalconst

"""
fid - city ID
da - xr.DataArray

Setting are
    clip = True (makes the read urban mask only covering the tiny area around the city instead of the whole US)
    opt = "tiff_3x"
    which = "core" (use this one because reviewer asking about within city plant types)
       ( = "both" means both the urban core and rural areas
         = "core" means only the urban area [Us]
         = "rural" means only the rural area )


For mask_low_evi_seasonal, maybe just set to winter         
"""

"""
Call chain logic

mask_water (internall reads NLCD data, which internally mask out pixels outside the city in my NLCD tiff files)
mask_crop (internall reads NLCD data, which internally mask out pixels outside the city in my NLCD tiff files)
mask_impervious (internally reads impervious data, which internally mask out pixels outside the city in my NLCD tiff files)

// skipping this step is okay; mask out few additional pixels
mask_low_evi_seasonal (internally reads EVI values, and for each season of the year, mask out where EVI < 0.05 in this season; can set season = "DJF" to only use winter EVI)
"""

"""
Apply these functions on the 30m data will make the consistent mask
"""


def get_mask(fid, which=["both", "core", "rural"], clip=[True, False], opt=["tiff", "tiff_3x"]):
    """ In this function, I get the urban mask geotiff to remove anything
        outside the urban boundary in NLCD & impervious area & EVI data.
        This function is used in read_nlcd() and read_impervious() below so I included it here.
    """

    if which == "core":
        ds = rio.open(os.path.join(path_intrim, "urban_mask", "US_urbanCluster.tif"))
        mask = ds.read()[0, :, :] == fid
        ds.close()
    else:
        if opt == "tiff":
            ds = rio.open(
                os.path.join(path_intrim, "urban_mask", "US_urbanCluster_merged.tif")
            )
            mask = ds.read()[0, :, :] == fid
            ds.close()
        else:
            ds = rio.open(
                os.path.join(path_intrim, "urban_mask", "city_boundary_3x_merged.tif")
            )
            mask = ds.read()[0, :, :] == fid
            ds.close()

        if which == "rural":
            ds = rio.open(
                os.path.join(path_intrim, "urban_mask", "US_urbanCluster.tif")
            )
            mask2 = ds.read()[0, :, :] == fid
            ds.close()
            mask = mask & (~mask2)

    if clip:
        # Clip mask to like GEE exports: one extra column to the right, one extra row to the bottom
        if opt == "tiff":
            ds = rio.open(
                os.path.join(path_intrim, "urban_mask", "US_urbanCluster_merged.tif")
            )
        else:
            ds = rio.open(
                os.path.join(path_intrim, "urban_mask", "city_boundary_3x_merged.tif")
            )
        temp = ds.read()[0, :, :] == fid
        retain = np.where(temp)
        ds.close()
        mask = mask[
               min(retain[0]): (max(retain[0]) + 2), min(retain[1]): (max(retain[1] + 2))
               ]

    return mask

def read_nlcd(fid, which=['core', 'rural', 'both'], opt=['tiff', 'tiff_3x']):
    """ Read NLCD data. These use the NLCD*.tif in the google drive folder I sent you. """
    path_intrim = "/Users/patrickchan/Downloads/"
    f = rio.open(os.path.join(path_intrim, 'gee_single', 'water and crop', "", f'NLCD_{fid:02d}.tif'))
    temp = f.read()
    f.close()

    dall = xr.DataArray(temp.reshape(8, -1, temp.shape[1], temp.shape[2]), dims=['year', 'band', 'row', 'col'],
                        coords={'year': [2001, 2004, 2006, 2008, 2011, 2013, 2016, 2019],
                                'band': [11, 12, 21, 22, 23, 24, 31, 41, 42, 43, 51, 52, 71, 72, 73, 74, 81, 82, 90,
                                         95],
                                # NLCD land cover class table; see https://developers.google.com/earth-engine/datasets/catalog/USGS_NLCD_RELEASES_2019_REL_NLCD#bands
                                'row': range(temp.shape[1]),
                                'col': range(temp.shape[2])})

    # remove the land cover types only in Alaska since they are irrelevant to CONUS
    dall = dall[:, ~np.isin(dall['band'].values, [51, 72, 73, 74]), :, :]

    # interpolate intermediate years
    dall = dall.interp({'year': range(2001, 2020)}, method='linear')

    # year 2020 = 2019
    dall = dall.interp({'year': range(2001, 2021)}, method='nearest', kwargs={'fill_value': 'extrapolate'})

    mask = get_mask(fid, which, True, opt)
    dall = dall.where(mask)
    return dall

def read_impervious(fid, which=["core", "rural", "both"], opt=["tiff", "tiff_3x"]):
    """ Read impervious area data. These use the impervious*.tif in the google drive folder I sent you. """
    f = rio.open(
        os.path.join(path_intrim, "gee_single", "Impervious", opt, f"impervious_{fid:02d}.tif"
        )
    )
    temp = f.read() / 100  # convert from percentage to fraction
    f.close()

    dall = xr.DataArray(
        temp,
        dims=["year", "row", "col"],
        coords={
            "year": [2001, 2004, 2006, 2008, 2011, 2013, 2016, 2019],
            "row": range(temp.shape[1]),
            "col": range(temp.shape[2]),
        },
    )

    # interpolate intermediate years
    dall = dall.interp({"year": range(2001, 2020)}, method="linear")

    # year 2020 = 2019
    dall = dall.interp(
        {"year": range(2001, 2021)},
        method="nearest",
        kwargs={"fill_value": "extrapolate"},
    )

    mask = get_mask(fid, which, True, opt)
    dall = dall.where(mask)
    return dall


####


####

def mask_water(fid, da, opt=["tiff", "tiff_3x"], which=["urban", "rural", "both"]):
    """Use NLCD to remove the pixels with > 40% water"""

    ###
    nlcd = read_nlcd(fid, which, opt).mean(
        axis=0
    )  # average percentage land cover over area over 2001-2020
    mask = nlcd.loc[11, :, :].values < 0.4
    if isinstance(da, xr.DataArray):
        da = da.where(mask)
    else:
        mask = np.where(mask, 1, np.nan)
        mask = np.broadcast_to(mask, da.shape)
        da = da * mask
    return da


def mask_impervious(
        fid, da, thres, opt=["tiff", "tiff_3x"], which=["urban", "rural", "both"]
):
    """Use NLCD to remove the pixels with impervious fraction > thres; thres is set to 0.8"""
    dall = read_impervious(fid, which, opt)
    mask = dall.mean(axis=0) < thres  # average impervious area over the time period
    if isinstance(da, xr.DataArray):
        da = da.where(mask)
    else:
        mask = np.where(mask, 1, np.nan)
        mask = np.broadcast_to(mask, da.shape)
        da = da * mask
    return da

def mask_crop(fid, da, opt=["tiff", "tiff_3x"], which=["urban", "rural", "both"]):
    """Use NLCD to remove the pixels with >50% crops """
    nlcd = read_nlcd(fid, which, opt).mean(
        axis=0
    )  # average percentage land cover over area over 2001-2020
    mask = nlcd.loc[82, :, :].values < 0.5
    if isinstance(da, xr.DataArray):
        da = da.where(mask)
    else:
        mask = np.where(mask, 1, np.nan)
        mask = np.broadcast_to(mask, da.shape)
        da = da * mask
    return da


# def mask_low_evi_seasonal(da_evi, fid, name="MOD09Q1G_EVI", extent=["tiff", "tiff_3x"], season=None):
#     """Use the MOD09Q1G_EVI_*_mask.tif files to remove the pixels that have seasonal mean EVI <= 0.05.
#     This is to ensure that the extreme event fell on enough vegetation
#     to have an impact on vegetation."""
#     h = rio.open(os.path.join(path_out, "veg", extent, f"{name}_{fid}_mask.tif"))
#     mask = h.read()[:4, :, :]
#     h.close()
#
#     if season is None:
#         for i, _ in enumerate(["DJF", "MAM", "JJA", "SON"]):
#             filt = da_evi["time"].to_index().to_period("Q-NOV").quarter == (i + 1)
#             da_evi[filt, :, :] = da_evi[filt, :, :].where(mask[i, :, :] > 0.05).values
#     else:
#         if season == "DJF":
#             quarter = 1
#         elif season == "MAM":
#             quarter = 2
#         elif season == "JJA":
#             quarter = 3
#         elif season == "SON":
#             quarter = 4
#         if isinstance(da_evi, xr.DataArray):
#             da_evi = da_evi.where(mask[quarter - 1, :, :] > 0.05)
#         else:
#             da_evi = np.where(mask[quarter - 1, :, :] > 0.05, da_evi, np.nan)
#     return da_evi

def write_final_mask(final_mask,NLCD_tif):
    """

    :param final_mask: it is a xarray
    :param NLCD_tif:
    :param out_tif
    :return:
    """
    #create tif

    out_tif = "NLCD_07_finalmask.tif" #the final output tif

    if os.path.exists(out_tif):
        os.remove(out_tif)
    else:
        pass

    NLCD_ds = gdal.Open(NLCD_tif, gdalconst.GA_ReadOnly)
    xsize = NLCD_ds.RasterXSize
    ysize = NLCD_ds.RasterYSize
    driver = gdal.GetDriverByName('GTiff')
    num_specific_classes = 19 ##should be consistent of years
    output_ds = driver.Create(out_tif, xsize, ysize,  num_specific_classes, gdalconst.GDT_Float32)

    #set projection and geotransform
    # Set the projection and geotransform for the output dataset
    geotransform = NLCD_ds.GetGeoTransform()
    output_ds.SetGeoTransform((geotransform[0], 1000, 0, geotransform[3], 0, -1000))
    output_ds.SetProjection(NLCD_ds.GetProjection())

    from tqdm import tqdm
    for year_index in tqdm(range(num_specific_classes), desc="Outer Loop"):
        print(f"processing the '{year_index+1}-th year'!")
        ##read the final mask
        mask_arr = final_mask[year_index,:,:].to_numpy()
        #the rows and columns of np array is inversed comared to that of gdal array
        mask_arr = np.reshape(mask_arr, (ysize, xsize))
        mask_arr
        ##write the final mask into tif
        output_ds.GetRasterBand(year_index + 1).WriteArray(mask_arr)

        # write band description
        output_ds.GetRasterBand(year_index + 1).SetDescription("mask_"+str(2000+year_index))

    print(f"'{out_tif}'is generated successfully!")

    output_ds = None
    NLCD_ds = None

    return True

#mask the final GLFC30_1k
def mask_final_GLFC_1k(final_mask, GLFC_1k,lc_list):
    """

    :param mask: the final mask --0,1
    :param GLFC_1k: the GLFC data product at 1km resolution and NA_projection_1k_2sp
                    for exmaple
    :return:
    """
    masked_GLFC_1k = "masked_NLCD_07_Boston.tif"
    if os.path.exists(masked_GLFC_1k):
        os.remove(masked_GLFC_1k)
    else:
        pass

    from osgeo import gdal, gdalconst

    ####open GLFC_1k get the data that will be masked
    GLFC_1k_ds = gdal.Open(GLFC_1k)
    #GetRasterBand(>=1),#get the first band,GetRasterBand(=1)

    if os.path.exists(masked_GLFC_1k):
        os.remove(masked_GLFC_1k)
    else:
        pass

    driver = gdal.GetDriverByName('GTiff')
    #input_ds.RasterXsize is the columns
    #input_ds.RasterYsize is the rows
    bands = GLFC_1k_ds.RasterCount
    output_ds = driver.Create(masked_GLFC_1k, GLFC_1k_ds.RasterXSize, GLFC_1k_ds.RasterYSize, bands, gdalconst.GDT_Float32)
    #set the number of bands as 1

    #when write an numpy into the "output_ds", its shape[0] should be "input_ds.RasterYsize"

    # Set the projection and geotransform for the output dataset
    geotransform = GLFC_1k_ds.GetGeoTransform()
    output_ds.SetGeoTransform((geotransform[0], 1000, 0, geotransform[3], 0, -1000))
    output_ds.SetProjection(GLFC_1k_ds.GetProjection())

    for year_index in range(1): #years =19, "1" here represents year 2000
        print(f"processing'{2000+year_index}-th year'")
        mask_arr = final_mask[year_index, :, :].to_numpy()

        # read each band of GLFC30_1k
        # each band represent percentage of a specified land cover
        for output_lc_index in range(bands): ## bands represents the 36 classes in a year
            print(f"processing'{output_lc_index+1}-th band'")
            # when mask=1, glfc pixel value will not be changed
            # when mask=1, glfc pixel value will change as nan
            GLFC_1k_ds
            GLFC_1k_arr = GLFC_1k_ds.GetRasterBand(output_lc_index+1).ReadAsArray(0, 0, GLFC_1k_ds.RasterXSize, GLFC_1k_ds.RasterYSize)
            GLFC_1k_arr
            band_desc = GLFC_1k_ds.GetRasterBand(output_lc_index+1).GetDescription()
            band_desc
            masked_glfc_arr = np.where(mask_arr==1,GLFC_1k_arr,np.nan)
            masked_glfc_arr
            output_ds.GetRasterBand(output_lc_index + 1).WriteArray(masked_glfc_arr)

            # write band description
            output_ds.GetRasterBand(output_lc_index + 1).SetDescription("masked_"+str(2000+year_index)+"_"+str(lc_list[output_lc_index]))

        print("finished")

    output_ds = None
    GLFC_1k_ds = None

    return True

if __name__ == '__main__':

    fid_list = [i for i in range(85)] #city id
    # fid = 1
    # path_intrim = "/Users/patrickchan/Downloads/"
    # nlcd = rio.open(os.path.join(path_intrim, 'gee_single', 'water and crop', "", f'NLCD_{fid:02d}.tif'))
    #
    # nlcd
    #
    # da = xr.DataArray(np.ones([19, nlcd.shape[0], nlcd.shape[1]]), dims=['year', 'row', 'col'],
    #                   coords={'year': range(2001, 2020), 'row': range(nlcd.shape[0]), 'col': range(nlcd.shape[1])})
    #
    # data = mask_water(fid, da, opt = "tiff_3x", which="urban")
    # data
    # print(data[0,:,:])
    # mask_mask_crop_m1 = data[0,:,:]
    # mask_mask_crop_m1
    # plt.imshow(data[0,:,:])

    fid = 7 #fid=7 Boston Area
    path_intrim = "/Users/patrickchan/Downloads/"
    nlcd = rio.open(os.path.join(path_intrim, 'gee_single', 'water and crop', "", f'NLCD_{fid:02d}.tif'))
    nlcd.shape[0]
    nlcd.shape[1]
    nlcd
    #2001-2019: 19 years
    da = xr.DataArray(np.ones([19, nlcd.shape[0], nlcd.shape[1]]), dims=['year', 'row', 'col'],
                      coords={'year': range(2001, 2020), 'row': range(nlcd.shape[0]), 'col': range(nlcd.shape[1])})

    mask_water_m1 = mask_water(fid, da, opt = "tiff_3x", which="urban")

    mask_crop_m1 = mask_crop(fid, da, "tiff_3x", which="urban")
    print("mask_crop \n")
    print(mask_crop_m1)

    mask_impervious_m1 = mask_impervious(fid, da, 0.8,  "", "urban")
    print("mask_impervious_m1 \n")
    print(mask_impervious_m1)

    mask_impervious_m1.values
    print("mask_impervious_m1 \n")
    print(mask_impervious_m1)

    # Set values to 1 where all three DataArrays have 1, keep NaN otherwise
    final_result = xr.where((mask_water_m1 == 1) & (mask_crop_m1 == 1) & (mask_impervious_m1 == 1), 1, np.nan)
    final_result
    print("logical results\n")
    print(final_result)

    #get tha mask and write it into a tif
    # NLCD_tif = os.path.join(path_intrim, 'gee_single', 'water and crop', "", f'NLCD_{fid:02d}.tif')
    # write_final_mask(final_result, NLCD_tif)

    GLFC_1k = "output_lc_perc_Project_Clip.tif"
    lc_cover_ls = [10, 11, 12, 20, 51, 52, 61, 62, 71, 72, 81, 82, 91, 92, 120, 121, 122, 130,
              140, 150, 152, 153, 181, 182, 183, 184, 185, 186, 187, 190, 200, 201, 202, 210, 0, 250]

    mask_final_GLFC_1k(final_result, GLFC_1k,lc_cover_ls)


