from datetime import datetime, timedelta
from netCDF4 import num2date, date2num
import numpy.ma as ma
from netCDF4 import Dataset
import argparse
import re
import numpy

# IMOS file format conversion to OceanSITES format
# Pete Jansen 2018-10-09

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

if not ncTime:
    print("time variable not found")
    # exit(-1)
    ncTime = {}
    ncTime[0] = nc.variables["TIME"]

time_deployment_start = nc.time_deployment_start
time_deployment_end = nc.time_deployment_end

tStart = parse(time_deployment_start, ignoretz=True)
tEnd = parse(time_deployment_end, ignoretz=True)

tStartnum = date2num(tStart, units=ncTime[0].units, calendar=ncTime[0].calendar)
tEndnum = date2num(tEnd, units=ncTime[0].units, calendar=ncTime[0].calendar)

maTime = ma.array(ncTime[0][:])
msk = (maTime < tStartnum) | (maTime > tEndnum)
maTime.mask = msk

dates = num2date(maTime.compressed(), units=ncTime[0].units, calendar=ncTime[0].calendar)

nc.close()


# create a new filename
# from:
# IMOS_<Facility-Code>_<Data-Code>_<Start-date>_<Platform-Code>_FV<File-Version>_ <Product-Type>_END-<End-date>_C-<Creation_date>_<PARTX>.nc
# to:
# OS_[PlatformCode]_[DeploymentCode]_[DataMode]_[PARTX].nc

# TODO: what to do with <Data-Code> with a reduced number of variables

splitPath = path_file.split("/")
fileName = splitPath[-1]
splitParts = fileName.split("_") # get the last path item (the file nanme), split by _

tStartMaksed = dates[0]
tEndMaksed = dates[-1]

fileProductTypeSplit = splitParts[6].split("-")
fileProductType = fileProductTypeSplit[0]

# could use the global attribute site_code for the product type

fileTimeFormat = "%Y%m%d"
ncTimeFormat = "%Y-%m-%dT%H:%M:%SZ"
sensor = fileProductTypeSplit[3]

# Generate OceanSITES file name
outputName = "OS" \
             + "_" + "SOTS" \
             + "_" + fileProductType + "-" + fileProductTypeSplit[1] \
             + "_D" \
             + "_" + sensor + "-" + splitParts[-1] \
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
                            'date_created', 'Mooring', 'project',
                            'deployment_code', 'Metadata_Conventions',
                            'instrument', 'disclaimer', 'distribution_statement',
                            'instrument_nominal_depth', 'standard_name_vocabulary',
                            'instrument_sample_interval',
                            'instrument_serial_number', 'institution', 'institution_address',
                            'quality_control_log', 'citation', 'contributor_role', 'data_centre', 'data_centre_email',
                            'history', 'acknowledgement', 'abstract', 'author', 'author_email', 'file_version'
                            'references']

# these attributes are actually part of the OceanSites list, so I took them off the BlackList
#, 'netcdf_version','geospatial_lat_max', 'geospatial_lat_min', 'geospatial_lon_max', 'geospatial_lon_min',
# 'geospatial_vertical_max', 'geospatial_vertical_min',]

# global attributes

dsIn = Dataset(path_file, mode='r')
for a in dsIn.ncattrs():
    if not (a in globalAttributeBlackList):
        print("Attribute %s value %s" % (a, dsIn.getncattr(a)))
        ncOut.setncattr(a, dsIn.getncattr(a))

# copy dimensions

for d in dsIn.dimensions:
    size = dsIn.dimensions[d].size
    if d == 'TIME':
        size = dates.shape[0]
    print("dimension", d, " shape ", size)
    ncOut.createDimension(d, size)

ncOut.createDimension('LATITUDE', 1)
ncOut.createDimension('LONGITUDE', 1)

history = dsIn.getncattr("history")
instrumentName = dsIn.getncattr("instrument")
serialNumber = dsIn.getncattr("instrument_serial_number")
deployment_start = dsIn.getncattr("time_deployment_start")
deployment_end = dsIn.getncattr("time_deployment_end")
deployment_voyage = dsIn.getncattr("voyage_deployment")
recovery_voyage = dsIn.getncattr("voyage_recovery")
deployment_splitparts = deployment_voyage.split("=")
recovery_splitparts = recovery_voyage.split("=")
platform_deployment_cruise_name = deployment_splitparts[-1]
platform_recovery_cruise_name = recovery_splitparts[-1]

if platform_deployment_cruise_name[:2]  == "IN":
    deployment_ship = "RV Investigator"
elif platform_deployment_cruise_name[:2]  == "SS":
    deployment_ship = "RV Southern Surveyor"

if platform_recovery_cruise_name[:2]  == "IN":
    recovery_ship = "RV Investigator"
elif platform_recovery_cruise_name[:2]  == "SS":
    recovery_ship = "RV Southern Surveyor"



dsIn.close()

ncOut.setncattr("time_coverage_start", dates[0].strftime(ncTimeFormat))
ncOut.setncattr("time_coverage_end", dates[-1].strftime(ncTimeFormat))
ncOut.setncattr("date_created", datetime.utcnow().strftime(ncTimeFormat))
ncOut.setncattr("history", history + '\n' + datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC : Convert from IMOS file : ") + fileName)
ncOut.setncattr("acknowledgement", "We acknowledge support from the following agencies: the Australian Antarctic Program Partnership (AAPP), the Antarctic Climate and Ecosystems Cooperative Research Centre (ACE CRC), the Integrated Marine Observing System (www.imos.org.au), University of Tasmania (UTAS), Bureau of Meteorology (BoM), the Marine National Facility (MNF) and the Australian Antarctic Division (AAD). We also acknowledge the support of the CSIRO Moored Sensor Systems team.")
ncOut.setncattr("data_type", "OceanSITES time-series data")
ncOut.setncattr("format_version", "1.3")
ncOut.setncattr("network", "IMOS")
# TODO: is this the right theme to use?
ncOut.setncattr("theme", "Global Ocean Watch")
ncOut.setncattr("summary", "Particle flux data from the Southern Ocean Time Series observatory in the Southern Ocean southwest of Tasmania.")
ncOut.setncattr("id", outputName)
ncOut.setncattr("sea_area", "Pacific Ocean")
ncOut.setncattr("array", "SOTS")
ncOut.setncattr("update_interval", "void")
ncOut.setncattr("creator_institution", "Commonwealth Scientific and Industrial Research Organisation (CSIRO)")
ncOut.setncattr("source", "subsurface mooring")
ncOut.setncattr("naming_authority", "OceanSITES")
ncOut.setncattr("data_assembly_center", "Australian Ocean Data Network (AODN)")
ncOut.setncattr("citation", "Any users of IMOS data are required to clearly acknowledge the source of the material derived from IMOS in the format: "
                            "\"Data was sourced from the Australian Integrated Marine Observing System (IMOS). IMOS is enabled by the National Collaborative Research Infrastructure Strategy (NCRIS). "
                            "It is operated by a consortium of institutions as an unincorporated joint venture, with the University of Tasmania as Lead Agent.\" "
                            "as well as \"PI Elizabeth Shadwick SOTS - SAZ. These data were collected and made freely "
                            "available by the OceanSITES program and the national programs that contribute to it. "
                            "[year-of-data-download], [Title], [Data access URL], accessed [date- of-access].\"")
ncOut.setncattr("Conventions", "CF-1.6 OceanSITES-1.3 NCADD-1.2.1")
ncOut.setncattr("license", "Follows CLIVAR (Climate Variability and Predictability) standards,cf. http://www.clivar.org/data/data_policy.php Data available free of charge. User assumes all risk for use of  data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data.")
ncOut.setncattr("contributor_name", "CSIRO; IMOS; ACE-CRC; MNF; AAPP")
ncOut.setncattr("processing_level", "data manually reviewed")
# TODO: can we say "excellent" for the SAZ data since it is QCed?
ncOut.setncattr("QC_indicator", "excellent")
ncOut.setncattr("instrument", instrumentName + "-" + serialNumber)
ncOut.setncattr("cdm_data_type", "Station")
ncOut.setncattr("platform_deployment_date", deployment_start)
ncOut.setncattr("platform_recovery_date", deployment_end)
ncOut.setncattr("creator_name", "Cathryn Wynn-Edwards")
ncOut.setncattr("creator_email", "cathryn.wynn-edwards@csiro.au")
ncOut.setncattr("creator_type", "person")
ncOut.setncattr("publisher_name", "Peter Jansen")
ncOut.setncattr("publisher_email", "peter.jansen@csiro.au")
ncOut.setncattr("netcdf_version", "netcdf-4 classic")
ncOut.setncattr("keywords_vocabulary", "GCMD Science Keywords")
ncOut.setncattr("keywords", "PARTICLE FLUX, CARBON, SEDIMENT COMPOSITION, INORGANIC CARBON,"
                            "MARINE GEOCHEMISTRY, ORGANIC CARBON, ORGANIC MATTER, SILICATE, CARBONATE")
ncOut.setncattr("platform_deployment_ship_name", deployment_ship)
ncOut.setncattr("platform_deployment_cruise_name", platform_deployment_cruise_name)
ncOut.setncattr("platform_recovery_ship_name", recovery_ship)
ncOut.setncattr("platform_recovery_cruise_name", platform_recovery_cruise_name)
ncOut.setncattr("references", "http://www.imos.org.au, data QC procedure document: http://dx.doi.org/10.26198/5dfad21358a8d, http://www.oceansites.org/")

# TODO: sort the attributes alphabetcially?



# copyData
#

# copy variable data from input files into output file

nc1 = Dataset(path_file, mode="r")

varList = nc1.variables

for v in varList:
    print("%s file %s" % (v, path_file))

    maVariable = nc1.variables[v][:]  # get the data
    varDims = varList[v].dimensions
    var_out = maVariable
    if 'TIME' in varDims:
        print("its a time variable shape ", var_out.shape, "dims", varDims, "len shape", len(var_out.shape))
        if varDims[0] != 'TIME':
            var_out = maVariable[:, ~msk]
            print("mask time dimension")
        elif len(var_out.shape) > 2:
            var_out = maVariable[~msk, :]
        else:
            var_out = maVariable[~msk]
    else:
        var_out = maVariable

    print("var out shape ", var_out.shape)

    #print(varDims)
    #print(maVariable.compressed().shape)

    # rename the _quality_control variables _QC
    varnameOut = re.sub("_quality_control", "_QC", v)

    fill = None
    try:
        fill = varList[v].fill_value
    except:
        pass
    ncVariableOut = ncOut.createVariable(varnameOut, varList[v].dtype, varDims, fill_value=fill)
    print("netCDF variable out shape", ncVariableOut.shape, "dims", varDims)
    # copy the variable attributes
    # this is ends up as the super set of all files
    for a in varList[v].ncattrs():
        if a != '_FillValue':
            print("%s Attribute %s = %s" % (varnameOut, a, varList[v].getncattr(a)))
            attValue = varList[v].getncattr(a)

            # Could restrict this just to the ancillary_variables attribute
            if isinstance(attValue, str):
                attValue = re.sub("_quality_control", "_QC", varList[v].getncattr(a))
            ncVariableOut.setncattr(a, attValue)

    ncVariableOut[:] = var_out

    # update the history attribute
    try:
        hist = nc1.history + "\n"
    except AttributeError:
        hist = ""

    ncOut.setncattr('history', hist + datetime.utcnow().strftime("%Y-%m-%d") + " : converted to oceanSITES format from file " + path_file)

nc1.close()

ncOut.close()

