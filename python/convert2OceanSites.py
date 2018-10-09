from datetime import datetime, timedelta
from netCDF4 import num2date, date2num
from netCDF4 import stringtochar
import numpy.ma as ma
import sys
from netCDF4 import Dataset
import numpy
import argparse

# IMOS file format convertion to OceanSITES format
# Pete Jansen 2019-10-09

from dateutil.parser import parse

parser = argparse.ArgumentParser()
parser.add_argument('file', help='input file name')
args = parser.parse_args()

path_file = args.file

# split this into   createCatalog - copy needed information into structure
#                   createTimeArray (1D, OBS) - from list of structures
#                   createNewFile
#                   copyAttributes
#                   updateAttributes
#                   copyData

#
# createCatalog - copy needed information into structure
#

print("input file %s" % path_file)

nc = Dataset(path_file, mode="r")

ncTime = nc.get_variables_by_attributes(standard_name='time')

time_deployment_start = nc.time_deployment_start
time_deployment_end = nc.time_deployment_end

tStart = parse(time_deployment_start)
tEnd = parse(time_deployment_end)

tStartnum = date2num(tStart.replace(tzinfo=None), units=ncTime[0].units, calendar=ncTime[0].calendar)
tEndnum = date2num(tEnd.replace(tzinfo=None), units=ncTime[0].units, calendar=ncTime[0].calendar)

maTime = ma.array(ncTime[0][:])
msk = (maTime < tStartnum) | (maTime > tEndnum)
maTime.mask = msk

dates = num2date(maTime, units=ncTime[0].units, calendar=ncTime[0].calendar)

nc.close()


# create a new filename
# from:
# IMOS_<Facility-Code>_<Data-Code>_<Start-date>_<Platform-Code>_FV<File-Version>_ <Product-Type>_END-<End-date>_C-<Creation_date>_<PARTX>.nc
# to:
# OS_[PlatformCode]_[DeploymentCode]_[DataMode]_[PARTX].nc

# TODO: what to do with <Data-Code> with a reduced number of variables

splitPath = path_file.split("/")
splitParts = splitPath[-1].split("_") # get the last path item (the file nanme), split by _

tStartMaksed = dates[0]
tEndMaksed = dates[-1]

fileProductTypeSplit = splitParts[6].split("-")
fileProductType = fileProductTypeSplit[0]

# could use the global attribute site_code for the product type

fileTimeFormat = "%Y%m%d"
ncTimeFormat = "%Y-%m-%dT%H:%M:%SZ"

outputName = "OS" \
             + "_" + "IMOS-EAC" \
             + "_" + fileProductType + "-" + fileProductTypeSplit[1] \
             + "_D" \
             + "_" + fileProductTypeSplit[2] + "-" + fileProductTypeSplit[3] + "m" \
             + ".nc"

print("output file : %s" % outputName)

ncOut = Dataset(outputName, 'w', format='NETCDF4')

#
# copyAttributes
#

# some of these need re-creating from the combined source data
globalAttributeBlackList = ['time_coverage_end', 'time_coverage_start',
                            'time_deployment_end', 'time_deployment_start',
                            'compliance_checks_passed', 'compliance_checker_version', 'compliance_checker_imos_version',
                            'date_created',
                            'deployment_code',
                            'geospatial_lat_max',
                            'geospatial_lat_min',
                            'geospatial_lon_max',
                            'geospatial_lon_min',
                            'geospatial_vertical_max',
                            'geospatial_vertical_min',
                            'instrument',
                            'instrument_nominal_depth',
                            'instrument_sample_interval',
                            'instrument_serial_number',
                            'quality_control_log',
                            'history', 'netcdf_version']


# global attributes

dsIn = Dataset(path_file, mode='r')
for a in dsIn.ncattrs():
    if not (a in globalAttributeBlackList):
        print("Attribute %s value %s" % (a, dsIn.getncattr(a)))
        ncOut.setncattr(a, dsIn.getncattr(a))

# copy dimensions

for d in dsIn.dimensions:
    ncOut.createDimension(d, dsIn.dimensions[d].size)

dsIn.close()

ncOut.setncattr("time_coverage_start", dates[0].strftime(ncTimeFormat))
ncOut.setncattr("time_coverage_end", dates[-1].strftime(ncTimeFormat))
ncOut.setncattr("date_created", datetime.utcnow().strftime(ncTimeFormat))
ncOut.setncattr("history", datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC : Convert from IMOS file"))

#
# copyData
#

# copy variable data from input files into output file

nc1 = Dataset(path_file, mode="r")

varList = nc1.variables

for v in varList:
    print("%s file %s" % (v, path_file))

    maVariable = nc1.variables[v][:]  # get the data
    varDims = varList[v].dimensions

    ncVariableOut = ncOut.createVariable(v, varList[v].dtype, varDims)
    # copy the variable attributes
    # this is ends up as the super set of all files
    for a in varList[v].ncattrs():
        print("%s Attribute %s value %s" % (v, a, varList[v].getncattr(a)))
        ncVariableOut.setncattr(a, varList[v].getncattr(a))

    ncVariableOut[:] = maVariable

nc1.close()

ncOut.close()

