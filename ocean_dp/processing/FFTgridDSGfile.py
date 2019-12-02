#!/usr/bin/python3

# readDSGfile
# Copyright (C) 2019 Peter Jansen
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


from netCDF4 import Dataset, num2date, chartostring
from dateutil import parser
from datetime import datetime
from datetime import timedelta

import numpy as np
import matplotlib.pyplot as plt

import finufftpy

import sys


def FFTgridDSGfile(netCDFfiles):
    ds = Dataset(netCDFfiles[1], 'a')

    vs = ds.get_variables_by_attributes(standard_name='sea_water_pressure_due_to_sea_water')
    vs = ds.get_variables_by_attributes(long_name='actual depth')

    pres_var = vs[0]
    pres = pres_var[:]

    #plt.plot(pres)
    #plt.show()

    temp_var = ds.variables["TEMP"]

    print("Read and convert time")
    time_var = ds.variables["TIME"]
    #time = num2date(time_var[:], units=time_var.units, calendar=time_var.calendar)
    #first_hour = time[0].replace(minute=0, second=0, microsecond=0)

    t = time_var[:] * 24
    hours = t

    # scale time to -pi to pi
    hours_min = np.min(hours)
    hours_max = np.max(hours)
    mid_hours = (hours_max + hours_min)/2

    print("time min, max", hours_min, hours_max, num2date(hours_max/24, units=time_var.units, calendar=time_var.calendar), num2date(hours_min/24, units=time_var.units, calendar=time_var.calendar), num2date(mid_hours/24, units=time_var.units, calendar=time_var.calendar))

    t2pi = 2*np.pi * (hours - mid_hours) / (hours_max - hours_min)

    # scale pressure to -pi to pi

    pres_min = np.min(pres)
    pres_max = np.max(pres)
    pres_mid = (pres_max + pres_min)/2

    print("pres min, max", pres_min, pres_max, pres_mid)

    d2pi = 2*np.pi * (pres - pres_mid) / (pres_max - pres_min)

    #plt.plot(t2pi)
    #plt.show()

    nt_points = int(hours_max - hours_min)
    nd_points = 20
    print(nt_points, nd_points)

    fft = np.zeros([nt_points, nd_points], dtype=np.complex128)

    print("Calc FFT")

    x = t2pi
    y = d2pi
    c = np.array(temp_var[:], dtype=np.complex128)
    c.real = temp_var[:]
    c.imag = 0
    print(c)

    print("size ", len(x), len(y), len(c))
    # finufftpy.nufft2d2(x, y, c, isign, eps, f, debug=0, spread_debug=0, spread_sort=2, fftw=0, modeord=0, chkbnds=1, upsampfac=2.0)
    res = finufftpy.nufft2d2(x, y, c, 0, 1e-12, fft)

    print("fft res", res)
    print(fft)

    for x in range(0, nd_points):
        plt.plot(abs(fft[:,x]))
    plt.grid(True)
    plt.show()

    index_var = ds.variables["instrument_index"]
    idx = index_var[:]
    instrument_id_var = ds.variables["instrument_id"]
    #time = num2date(time_var[:], units=time_var.units, calendar=time_var.calendar)

    #print(idx)
    i = 0
    for x in chartostring(instrument_id_var[:]):
        #print (i, x, time[idx == 1], pres[idx == i])
        #plt.plot(time[idx == i], pres[idx == i])  # , marker='.'
        i += 1

    plt.gca().invert_yaxis()
    plt.grid(True)

    # close the netCDF file
    ds.close()

    plt.show()


if __name__ == "__main__":
    FFTgridDSGfile(sys.argv)
