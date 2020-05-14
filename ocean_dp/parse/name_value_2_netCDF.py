#!/usr/bin/python3

# raw2netCDF
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

import sys
import re
import os

from datetime import datetime
from netCDF4 import date2num, num2date
from netCDF4 import Dataset
import numpy as np

# parsers need to output
#  instrument
#  instrument_serial_number
#  time_coverage_start
#  time_coverage_end
# optional
#  date_created
#  history
#
# convert time to netCDF cf-timeformat (double days since 1950-01-01 00:00:00 UTC)

#
# parse the file
#

# 2010-03-05 09:15:00,BV=14.29650497,ZACCEL=-8.105728149,PAR=0,LOAD=1.5
# 2010-03-05 09:45:00,CHL=0,PAR=6
# 2010-03-05 10:15:00,BV=14.29650497,ZACCEL=-8.10522747,PAR=0,LOAD=3
# 2010-03-05 10:45:00,CHL=0,PAR=9

# 2012-07-06 20:12:42 INFO: 7 done time 96 ,BV=NAN ,PT=NAN ,OBP=31.24 ,OT=9.24 ,CHL=106 ,NTU=174 ,PAR=0.4665374 ,meanAccel=-9.764508 ,meanLoad=NAN
# 2012-07-06 21:12:54 INFO: 7 done time 109 ,BV=NAN ,PT=NAN ,OBP=31.24 ,OT=9.33 ,CHL=231 ,NTU=4136 ,PAR=0.4664596 ,meanAccel=-9.767615 ,meanLoad=NAN
# 2012-07-06 22:12:57 INFO: 7 done time 112 ,BV=NAN ,PT=NAN ,OBP=31.23 ,OT=9.4 ,CHL=325 ,NTU=3957 ,PAR=0.4665374 ,meanAccel=-9.770709 ,meanLoad=NAN


def parse(file, name):

    hdr = True
    dataLine = 0
    number_samples_read = 0
    nVars = 0
    data = []
    raw = []
    ts = []
    settings = []
    instrument_model = 'LI-190SA'
    instrument_serial_number = 'unknown'
    sep = ','

    filepath = file[0]

    with open(filepath, 'r', errors='ignore') as fp:
        line = fp.readline()
        cnt = 1
        while line:
            #print("Line ", line)
            lineSplit = line.strip().split(sep)
            #print('Split ', lineSplit)

            t = datetime.strptime(lineSplit[0][0:19], '%Y-%m-%d %H:%M:%S')

            # turn name=value into dict, must be a better way
            names = []
            values = []
            for i in lineSplit[1:]:
                nv = i.split("=")
                if len(nv) > 1:
                    names.append(nv[0].strip(" "))
                    values.append(nv[1].strip(" "))
            name_value = dict(zip(names, values))
            #print(t, name_value)

            if name in names:
                try:
                    data.append(float(name_value[name]))
                    ts.append(t)
                    number_samples_read = number_samples_read + 1

                except ValueError:
                    pass

            line = fp.readline()

    # trim data
    print("samplesRead %d" % (number_samples_read))

    if number_samples_read == 0:
        return

    print('instrument ', instrument_model, 'serial', instrument_serial_number)

    #
    # build the netCDF file
    #

    ncTimeFormat = "%Y-%m-%dT%H:%M:%SZ"

    outputName = filepath + ".nc"

    print("output file : %s" % outputName)

    ncOut = Dataset(outputName, 'w', format='NETCDF4')

    ncOut.instrument = 'LI-COR ; ' + instrument_model
    ncOut.instrument_model = instrument_model
    ncOut.instrument_serial_number = instrument_serial_number

    for s in settings:
        ncOut.setncattr("comment_file_settings_" + s[0], s[1])

    #     TIME:axis = "T";
    #     TIME:calendar = "gregorian";
    #     TIME:long_name = "time";
    #     TIME:units = "days since 1950-01-01 00:00:00 UTC";

    tDim = ncOut.createDimension("TIME", number_samples_read)
    ncTimesOut = ncOut.createVariable("TIME", "d", ("TIME",), zlib=True)
    ncTimesOut.long_name = "time"
    ncTimesOut.units = "days since 1950-01-01 00:00:00 UTC"
    ncTimesOut.calendar = "gregorian"
    ncTimesOut.axis = "T"
    ncTimesOut[:] = date2num(ts, units=ncTimesOut.units, calendar=ncTimesOut.calendar)

    ncVarOut = ncOut.createVariable(name, "f4", ("TIME",), fill_value=np.nan, zlib=True) # fill_value=nan otherwise defaults to max
    ncVarOut.long_name = 'PAR sensor Voltage'
    ncVarOut.units = 'mV'
    ncVarOut.sensor_SeaVoX_L22_code = 'SDN:L22::TOOL0193'
    ncVarOut.comment_sensor_type = 'cosine sensor'
    ncVarOut[:] = data

    #ncVarOut = ncOut.createVariable('PAR_RAW', "f4", ("TIME",), fill_value=np.nan, zlib=True) # fill_value=nan otherwise defaults to max
    #ncVarOut.units = '1'
    #ncVarOut[:] = raw

    ncOut.setncattr("time_coverage_start", num2date(ncTimesOut[0], units=ncTimesOut.units, calendar=ncTimesOut.calendar).strftime(ncTimeFormat))
    ncOut.setncattr("time_coverage_end", num2date(ncTimesOut[-1], units=ncTimesOut.units, calendar=ncTimesOut.calendar).strftime(ncTimeFormat))
    ncOut.setncattr("date_created", datetime.utcnow().strftime(ncTimeFormat))
    ncOut.setncattr("history", datetime.utcnow().strftime("%Y-%m-%d") + " created from file " + os.path.basename(filepath))

    ncOut.close()

    return outputName


if __name__ == "__main__":
    parse(sys.argv[1:], "PAR")

