# -*- coding: utf-8 -*-
"""
Calculate multidimentional time series.

@author: Michael

Provide necessary calculation methods for heatwave study.
"""

import netCDF4 as nc
import numpy as np
import xarray as xr
import scipy.stats as sps
import concurrent.futures
import os

def get_percentile(Filename, percentile):
    """
    

    Parameters
    ----------
    Filename : 3D NetCDF Dataset
        Data source containing a ['time'] dimension.
    percentile : Float, 0<value<1
        Desired percentile (left to right).

    Returns
    -------
    p : 2D NumPy array
        90th percentile of Filename values across time.

    """

    ds = xr.open_dataset(Filename)
    tmax = ds['tmax']
    
    #Getting necessary descriptives
    mu = tmax.mean(dim=['time'], skipna=True)
    sigma = tmax.std(dim=['time'], skipna=True)
    
    #Generate distribution
    dist = sps.norm(loc=mu, scale=sigma)
    p = dist.ppf(percentile)

    return p

def get_case(directory, p, mindays, isevent=False, display=False):
    """
    

    Parameters
    ----------
    directory : LIST
        List of COMPLETE directories.
    p : NUMPY ARRAY
        Percentile.
    mindays : INT
        Minimum day to be considered heat wave.

    Returns
    -------
    yearboollist : LIST of TUPLES
        List of tuples (Year, NumPy bool arrays).

    """
    yearboollist = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=(16)) as e:
        futures = [e.submit(caseof, dataset, p, mindays, isevent, display) for dataset in directory]
        for group in concurrent.futures.as_completed(futures):
            yearboollist.append(group.result())

    if isevent:
        print('Event count finished for list of', len(yearboollist), 'arrays.')
    else:
        print('Computation process finished for list of', len(yearboollist), 'arrays.')
    return sorted(yearboollist)

#====Supplementary methods====#

def caseof(path, percentileval, minday, casecount=True, display=False):
    """
    

    Parameters
    ----------
    path : String
        Must point to a 3D netCDF file following var[time,lat,lon].
    percentileval : NumPy Array
        Must be a 2D array.

    Returns
    -------
    casecount=False
    (int year,arraydims[time,lat,lon])
    NumPy Array with no Attributes
        Consisting of 3 dimensions and bool values.
    casecount=True
    (int year,arraydims[lat,lon])

    """

    print("Acquiring Heat Day")
    # file是个字典，其中"tmax"键值对是一个三维数组
    file = nc.Dataset(path)
    isheatday = np.empty((file.variables['time'].size, 360, 720))
    isheatday[:] = np.nan

    # i至少为360，j至少为720
    i = file.variables['lat'].size
    v = file.variables['lon'].size
    t = file.variables['time'].size

    data = file.variables['tmax'][:t, :i, :v]

    # index为提取所有not nan值的索引
    index = ~np.isnan(data[0])

    # 得到在percentileval的需要需要比较的对应值
    cmp_values_1d = percentileval[:i, :v][index]

    # 得到所有需要与percentileval比较的对应的数据
    cmp_values_2d = data[:, index]

    # 将2组数据进行比较，得到比较结果
    cmp_ret_2d = (cmp_values_2d > cmp_values_1d)

    # 进行赋值并返回
    isheatday[:, :i, :v][:, index] = cmp_ret_2d

    def __lessdayremover(origset, minday):

        time = 0
        thrownvalue = 0
        totaltrue = 0

        dataset = origset
        eventcount = np.zeros((360, 720))

        for lat in range(dataset[1, :, 1].size):
            for lon in range(dataset[1, 1, :].size):

                counter = 0
                for time in range(dataset[:, 1, 1].size):
                    """
                        Time pointer on each lat lon.
                        If is True, add 1 to counter;
                        If is False, check if counter<minday, if so, delete all true counted
                        Each heat case end with counter>3, then counter=0
                        Nan will be skipped
                        """
                    if dataset[time, lat, lon]:
                        totaltrue += 1
                        counter += 1
                    # False counter<minday
                    elif dataset[time, lat, lon] == False:
                        if counter < minday:
                            for i in range(1, counter + 1):
                                dataset[time - i, lat, lon] = False
                                totaltrue -= 1
                                thrownvalue += 1
                            counter = 0
                        # False counter>=minday, eventcount[lat lon]++
                        else:
                            eventcount[lat, lon] += 1
                            counter = 0
        if display:
            print("")
            print(path, ":")
            print("Values Thrown: ", thrownvalue, " out of ", dataset.size)
            print("Current True Values: ", totaltrue, " out of ", dataset.size)

        print()
        print("Array ready")

        res_mask = np.ma.getmask(data)

        if casecount:
            return np.ma.masked_array(eventcount, mask=res_mask[0])
        else:
            return np.ma.masked_array(dataset, mask=res_mask)

    print("Analyzing Data")

    file_name = os.path.basename(path)
    return (int(file_name.split(".")[1]), __lessdayremover((np.array(isheatday, dtype=np.bool)), minday))

