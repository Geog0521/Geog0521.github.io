#include <iostream>
#include <math.h>
#include "gdal_priv.h"
#include <omp.h>
#include <time.h>
#include <string>
#include <sstream>
#include <list>

// #include <algorithm>
// #include <dirent.h>
// #include <stdlib.h>
#include <string.h>
// #include <unistd.h>
// #include <vector>
#include <string>
using namespace std;
// using namespace cv;

int One_to_four(int &res, int &Mean, int *pixel_data_window, int nPixels)
{
    res = 0;
    Mean = 0;
    double Sum = 0;
    int Min = 999999999;
    int Max = -999999999;
    for (size_t i = 0; i < nPixels; i++)
    {
        Sum +=pixel_data_window[i];
        Min = Min < pixel_data_window[i] ? Min: pixel_data_window[i];
        Max = Max > pixel_data_window[i] ? Max: pixel_data_window[i];
    }
    Mean = round(Sum / nPixels);
    double Tempvalue = (Sum - Max - Min) / 2;
    int xa = floor(Tempvalue);
    int xb = ceil(Tempvalue);
    int da = xa - Min;
    int db = Max - xb;
    int d = min(da, db);
    if (d == 0)
    {
        if (da == db)
        {
            if (xa == xb)
                res = 1;
            else
                res = 6;
        }
        else
        {
            if (xa == xb)
                res = 4;
            else
                res = 12;
        }
    }
    else
    {
        if (da == db)
        {
            if (xa == xb)
                res = 12 + 24 * (d - 1) + 6;
            else
                res = 24 + 24 * (d - 1) + 6;
        }
        else
        {
            if (xa == xb)
                res = 12 + 24 * (d - 1) + 12;
            else
                res = 24 + 24 * (d - 1) + 12;
        }
    }
    return 0;
}

int Iteration(unsigned short *new_grid, double &res_value, unsigned short * banddata, int mybase, int rows_n, int cols_n)
{
    double *result = new double[(rows_n-1) * (cols_n-1)]();
    #pragma omp parallel for
    for (size_t rn = 0; rn < (rows_n-1); rn++)
    {
        // #pragma omp parallel for
        for (size_t cn = 0; cn < (cols_n-1); cn++)
        {
            int pixel_data_window[4] = {banddata[rn*cols_n + cn],banddata[(rn+1)*cols_n + cn],banddata[rn*cols_n + cn + 1],banddata[(rn+1)*cols_n + cn+1]};
            int temp = 0;
            int Mean = 0;
            One_to_four(temp, Mean, pixel_data_window, 4);
            new_grid[rn*(cols_n-1) + cn] = Mean;
            result[rn*(cols_n-1) + cn] = log((double)temp) / log((double)mybase);            
        }
    }
    for (size_t i = 0; i < (rows_n-1) * (cols_n-1); i++)
    {
        res_value += result[i];
    }
    delete[] result;
    result = NULL;
    return 0;
}

double GetBoltzmann_resampling(unsigned short *banddata, int mybase, bool relative, bool Normalization, int rows, int cols)
{
    double *res_value_list = new double[rows * cols]();
    int iCount = 0;
    unsigned short *pusTmpBanddata = new unsigned short[rows * cols]();
    memcpy(pusTmpBanddata, banddata, rows * cols * sizeof(unsigned short));
    int rows_n = rows;
    int cols_n = cols;
    unsigned short *new_grid = NULL;
    while (rows_n > 1 && cols_n > 1)
    {
        // cout << rows_n << endl;
        new_grid = new unsigned short[(rows_n - 1) * (cols_n - 1)]();
        double res_value = 0.0;
        Iteration(new_grid, res_value, pusTmpBanddata, mybase, rows_n, cols_n);
        res_value_list[iCount++] = res_value;
        if (relative)
            break;
        else
        {
            // cout<< rows_n<<endl;
            delete []pusTmpBanddata;
            pusTmpBanddata = new unsigned short[(rows_n - 1) * (cols_n - 1)]();
            memcpy(pusTmpBanddata, new_grid, (rows_n - 1) * (cols_n - 1) * sizeof(unsigned short));
            rows_n -= 1;
            cols_n -= 1;
            delete[] new_grid; 
        }
    }
    if (relative)
    {
        if (Normalization)
            return res_value_list[0]/(rows*cols);
        else
            return res_value_list[0];
    }
    else
    {
        double dSumRes = 0.0;
        for (size_t i = 0; i < iCount; i++)
        {
            dSumRes += res_value_list[i];
        }
        
        if (Normalization)
            return dSumRes/(rows*cols);
        else
            return dSumRes;
    }
    delete[]res_value_list;
    res_value_list = NULL;
    delete[]pusTmpBanddata;
    pusTmpBanddata = NULL;
    delete[]new_grid;
    new_grid = NULL;
}

int Iteration_Aggregation(unsigned short *new_grid2, double &RBE, unsigned short * banddata, int mybase, int rows_n, int cols_n)
{
    RBE = 0.0;
    for (size_t rn = 0; rn < rows_n/2; rn++)
    {
        for (size_t cn = 0; cn < cols_n/2; cn++)
        {
            int pixel_data_window[4] = {banddata[2 * rn * cols_n + 2 * cn], banddata[(2 * rn + 1) * cols_n + 2 * cn], banddata[2 * rn * cols_n + 2 * cn + 1], banddata[(2 * rn + 1) * cols_n + 2 * cn + 1]};
            int temp = 0;
            int Mean = 0;
            One_to_four(temp, Mean, pixel_data_window, 4);            
            RBE += log(temp) / log(mybase);
            new_grid2[rn * (cols_n / 2) + cn] = Mean;
        }
    }
    return 0;
}

double GetBoltzmann_Aggregation(unsigned short *banddata, int mybase, bool relative, bool Normalization, int rows, int cols)
{
    double ABE = 0.0;
    int nExp = ceil(min(log(rows)/log(2), log(cols)/log(2)));
    double *RBEs_list = new double[nExp]();
    int iExp = 0;
    unsigned short *pusTmpBanddata = new unsigned short[rows * cols]();
    memcpy(pusTmpBanddata, banddata, rows * cols * sizeof(unsigned short));
    int _nRows = rows;
    int _nCols = cols;

    unsigned short *new_grid2 = NULL;
    while (_nRows > 1 && _nCols > 1)
    {
        new_grid2 = new unsigned short[(_nRows / 2) * (_nCols / 2)]();
        double RBE = 0.0;
        Iteration_Aggregation(new_grid2,RBE,pusTmpBanddata,mybase,_nRows, _nCols);
        RBEs_list[iExp++] = RBE;
        if (relative)
            break;
        else
        {
            // cout<<(_nRows / 2) * (_nCols / 2)<<endl;
            delete []pusTmpBanddata;
            pusTmpBanddata = new unsigned short[(_nRows / 2) * (_nCols / 2)]();
            memcpy(pusTmpBanddata, new_grid2, (_nRows / 2) * (_nCols / 2) * sizeof(unsigned short));
            _nRows /= 2;
            _nCols /= 2;
            delete[] new_grid2;
        }
    }
    if (relative)
    {
        if (Normalization)
        {
            return RBEs_list[0] / (rows*cols);
        }
        else
        {
            return RBEs_list[0];
        }
    }
    else
    {
        for (size_t i = 0; i < iExp; i++)
            ABE += RBEs_list[i];
        if (Normalization)
        {
            return ABE / (rows*cols);
        }
        else
        {
            return ABE;
        }
    }
    
    // delete[]RBEs_list;
    // RBEs_list = NULL;
    // delete[]pusTmpBanddata;
    // pusTmpBanddata = NULL;
    // delete[]new_grid2;
    // new_grid2 = NULL;
}

int main(int argc, char *argv[])
{
    cout<<"Computing file "<<argv[1]<<endl;
    double dTStart = omp_get_wtime();
    double sum_entropy = 0;
    string str_band_lis;
    int mybase = atoi(argv[2]);  
    bool relative = atoi(argv[3]); 
    bool Normalization = atoi(argv[4]);
    bool all_bands_bool = atoi(argv[5]); // whether or not all of bands are used. O for No; 1 for yes
    str_band_lis = atoi(argv[6]); //bands_list
    //cout<<"Computing Boltzmann entropies of bands: "<<argv[6]<<endl;
    GDALAllRegister();
    CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
    GDALDataset* ImageMulti;
    int BandNum = 0;
    int Width = 0;
    int Height = 0;
    
    ImageMulti = (GDALDataset *)GDALOpen(argv[1], GA_ReadOnly);
    if (ImageMulti == nullptr)
        return false;
    BandNum = ImageMulti->GetRasterCount();
    Width = ImageMulti->GetRasterXSize();
    Height = ImageMulti->GetRasterYSize();
    unsigned short **InBands = new unsigned short*[BandNum]();
    //get all bands first and store them in a multi-dimensional array
    for (size_t iBand = 0; iBand < BandNum; iBand++)
    {
        InBands[iBand] = new unsigned short[Width * Height]();
        GDALRasterBand *pInBands = ImageMulti->GetRasterBand(iBand + 1);
        pInBands->RasterIO(GF_Read, 0, 0, Width, Height, InBands[iBand], Width, Height, GDT_UInt16, 0, 0);
    }
    if (all_bands_bool) //we calculate Boltzmann entropies of all bands
    {
        for (int iBand = 0; iBand < BandNum; iBand++)
        {
           printf("Running the %dnd band: ", iBand + 1);
           double tStartBand = omp_get_wtime();
        
           double dResult = GetBoltzmann_resampling(InBands[iBand], mybase, relative, Normalization, Height, Width);
           // double dResult = GetBoltzmann_Aggregation(InBands[iBand], mybase, relative, Normalization, Height, Width);
           sum_entropy = sum_entropy + dResult;
           printf("%f",dResult);
           cout<<" bits/pixel, time: "<<(omp_get_wtime() - tStartBand)<<" s"<<endl;
        }
    cout<<"Boltzmann entropy: "<<(sum_entropy)<<" bits/pixel"<<endl;
    cout<<"Total running time: "<<(omp_get_wtime() - dTStart)<<" s"<<endl;
    }
    else //we calculate Boltzmann entropies of specific bands
    {
      //read "list" str and convert it into a C++-based list
      cout<<"Computing Boltzmann entropies of bands: "<<argv[6]<<endl;
      list<int> lst;
      size_t pos = 0;
      string str_band_lis=argv[6];
      while ((pos = str_band_lis.find(" ")) != string::npos) {
        string token = str_band_lis.substr(0, pos);
        int num = stoi(token);
        lst.push_back(num);
        //cout<<"number: "<<num<<endl;
        str_band_lis.erase(0, pos + 1);
    }
    //cout<<"last str: "<<str_band_lis<<endl;
   
    //deal with the last number
    int num2 = stoi(str_band_lis);
    lst.push_back(num2);
     //for (int x : lst){
        //printf("Running the %d: ", x);
    //}
    for (int x : lst){
        printf("Running the %dnd band: ", x);
        
        if (x>BandNum)
        {
          cout<<"the input is larger than the total number of band!"<<endl;
          std::exit(0);
        }
        else
        {
          double tStartBand = omp_get_wtime();
          double dResult = GetBoltzmann_resampling(InBands[x-1], mybase, relative, Normalization, Height, Width);
          // double dResult = GetBoltzmann_Aggregation(InBands[iBand], mybase, relative, Normalization, Height, Width);
          sum_entropy = sum_entropy + dResult;
          printf("%f",dResult);
          cout<<" bits/pixel, time: "<<(omp_get_wtime() - tStartBand)<<" s"<<endl;
        }
      }
    }
    return 0;
}
