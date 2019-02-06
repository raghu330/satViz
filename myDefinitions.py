#!/usr/bin/env python
# coding=utf-8
__author__ = "raghav"

'''
What does this code do?
Container for all custom definitions.

HISTORY:
#1. 09-Jul-2018: Mihir version
'''

#- Start definitions
# definition-1
def dateTimeManipulation():
    """
    What does this definition do?
    Manipulates date time strings as per needs

    :return: dateSnap, dt, yr, mm, dd, tm, hrs, mins, dt_2, yr_2, mm_2, dd_2, tm_2, hrs_2, mins_2
    """
    # import necessary modules.
    import os, sys
    from datetime import datetime, timedelta

    # get current time in IST
    nowTm = str(datetime.now()).split()
    nowTime = nowTm[1].split(':')
    nowHr = nowTime[0]
    nextDay = str(datetime.now() + timedelta(days = 1)).split()
    nextDayDate = nextDay[0].split('-')
    nextDateSnap = nextDayDate[0] + nextDayDate[1] + nextDayDate[2]

    # Current System Date and Time manipulation in UTC
    dt = str(datetime.utcnow()).split()
    yr = dt[0].split('-')[0]
    mm = dt[0].split('-')[1]
    dd = dt[0].split('-')[2]
    tm = dt[1].split(':')
    hrs = tm[0]
    mins = tm[1]
    sec = tm[2]

    # Past 2 days test..
    dt_2 = str(datetime.utcnow() - timedelta(hours = 0.4)).split()
    yr_2 = dt_2[0].split('-')[0]
    mm_2 = dt_2[0].split('-')[1]
    dd_2 = dt_2[0].split('-')[2]
    tm_2 = dt_2[1].split(':')
    hrs_2 = tm_2[0]
    mins_2 = tm_2[1]
    sec_2 = tm_2[2]

    # Set basic string variables
    dateSnap = yr + mm + dd
    return dateSnap, dt, yr, mm, dd, tm, hrs, mins, dt_2, yr_2, mm_2, dd_2, tm_2, hrs_2, mins_2
# end definition

# definition-2
def mkFldrs(dateSnap):
    # """
    # What does this definition do?
    # This definition simply creates the necessary folders as per need.
    #
    # :param dateSnap: Date (YYYYMMDD) name as string
    # :return: basDir, datDir, outDir, logDir, webDir, geoTdir, exeDir, GSHHS_ROOT
    # """

    import os, ConfigParser

    # Set up basic data folders
    config = ConfigParser.ConfigParser()
    configFilePath = r'/home/satviz/Work/scripts/setFldrs.txt'
    config.readfp(open(configFilePath))

    # Parse folder paths from configuration file
    basDir = config.get('paths', 'basDir')
    datDir = config.get('paths', 'datDir')
    outDir = config.get('paths', 'outDir')
    logDir = config.get('paths', 'logDir')
    geoTdir = config.get('paths', 'geoTdir')
    webDir = config.get('paths', 'webDir')
    GSHHS_ROOT = config.get('paths', 'GSHHS_dir')
    srcDir = config.get('paths', 'srcDir')
    exeDir = config.get('paths', 'exeDir')
    tmpDir = config.get('paths', 'tmpDir')

    # Set-up the appropriate paths
    datDir = datDir + dateSnap + '/'
    outDir = outDir + dateSnap + '/'
    logDir = logDir + dateSnap + '/'
    geoTdir = geoTdir + dateSnap + '/'
    webDir = webDir + dateSnap + '/'
    msg1Src = srcDir + dateSnap + '/' + 'MSG1/H-000/'

    # Create directories and log files if not present.
    if not os.path.exists(datDir):
        os.makedirs(datDir)
        print("mkFldr() says: Created a new folder at %s " % datDir)
    else:
        print("datDir was already created! Check %s " % datDir)
    # end if-not condition

    if not os.path.exists(outDir):
        os.makedirs(outDir)
        print("mkFldr() says: Created a new folder at %s " % outDir)
    else:
        print("outDir was already created! Check %s " % outDir)
    # end if-not condition

    if not os.path.exists(logDir):
        os.makedirs(logDir)
        print("mkFldr() says: Created a new folder at %s " % logDir)
    else:
        print("logDir was already created! Check %s " % logDir)
    # end if-not condition

    if not os.path.exists(webDir):
        os.makedirs(webDir)
        print("mkFldr() says: Created a new folder at %s " % webDir)
    else:
        print("webDir was already created! Check %s " % webDir)
    # end if-not condition

    if not os.path.exists(geoTdir):
        os.makedirs(geoTdir)
        print("mkFldr() says: Created a new folder at %s " % geoTdir)
    else:
        print("geoTDir was already created! Check %s " % geoTdir)
    # end if-not condition

    # Create log files for finished & decompressed time slots
    finTmsFile = logDir + "finishedTimeSlots_" + dateSnap + ".txt"
    if not os.path.exists(finTmsFile):
        open(finTmsFile, 'w').close()
    # end if-not condition

    deCmprssList = logDir + 'decompressedTimeSlots_' + dateSnap + '.txt'
    if not os.path.exists(deCmprssList):
        open(deCmprssList, 'w').close()
    # end if-not condition

    return basDir, datDir, outDir, logDir, webDir, geoTdir, msg1Src, exeDir, GSHHS_ROOT, tmpDir
# end definition

# definition-3
def embellish(basDir, GSHHS_ROOT, imgStr, ii, dateSnap, timeSnap):
    """
    What does this definition do?
    Embellishes the image with custom graphics

    :param basDir: Base directory path
    :param GSHHS_ROOT: GSHHS installation folder
    :param imgStr: Complete path of the output image data as string
    :param ii: Channel name as string
    :param dateSnap: Date (YYYYMMDD) name as string
    :param timeSnap: Time (HHMM) name as string
    :return: Image with all the decorations..

    References:
    [1]. https://stackoverflow.com/questions/18522295/python-pil-change-greyscale-tif-to-rgb
    """

    import os, sys
    import aggdraw, PIL
    from PIL import Image, ImageFont
    from pydecorate import DecoratorAGG as dag
    from pycoast import ContourWriter
    from satpy.resample import get_area_def

    img = Image.open(imgStr, mode = "r")

    # Convert the image into RGB if it's not
    if ("L" or "A" in img.getbands()):
        rgbimg = Image.new("RGBA", img.size)
        rgbimg.paste(img)
        img = rgbimg
    else:
        print("\n It was already an RGB image, so no need to convert!")
    # end if-confition

    # Add logos
    dc = dag(img)
    # dc.add_logo(basDir + 'logo/NCMRWF.png', height = 75.0)
    dc.add_logo(basDir + 'logo/NCMRWF.png', height = 75.0, bg='yellow')

    # Add basic text information
    capStr = ii
    textStr = "MSG-1: SEVIRI [" + capStr + "]\n" + dateSnap + '\n' + timeSnap + ' UTC'
    # fontClr = aggdraw.Font((255, 255, 255), "/usr/share/fonts/truetype/DejaVuSerif.ttf", size = 18)
    fontClr = aggdraw.Font('black', "/usr/share/fonts/truetype/DejaVuSerif.ttf", size = 18)
    # fontClr = aggdraw.Font('blue', "/usr/share/fonts/liberation/LiberationMono-Bold.ttf", size = 18)
    dc.add_text(textStr, font = fontClr)

    # Adding partner's logos
    dc.align_right()
    dc.add_logo(basDir + 'logo/imdlogo.jpg', height=75.0)
    dc.align_left()
    dc.align_bottom()
    dc.add_logo(basDir + 'logo/eumetsatLogo.jpg', height=50.0)
    dc.align_right()
    dc.add_logo(basDir + 'logo/pytroll-logo.png', height = 50.0)

    # Add coastline, grid to the image
    # Set up projection info for coastlines
    proj4_str = '+proj=merc +lon_0=75.0 + lat_0=17.5 + lat_ts=17.5 +ellps=WGS84'
    area_extent = (-3339584.72, -1111475.10, 3339584.72, 5591295.92)
    area_def = (proj4_str, area_extent)

    # Add the shape file directory or the GSHHS root path
    cw = ContourWriter(GSHHS_ROOT)
    cw.add_coastlines(img, area_def, resolution = 'h', level = 1)
    cw.add_shapefile_shapes(img, area_def, GSHHS_ROOT + 'India/Admin2.shp', outline = 'white')

    # Add gridlines
    fontClr2 = ImageFont.truetype("/usr/share/fonts/truetype/DejaVuSerif.ttf", 14)
    # fontClr2 = ImageFont.truetype("/usr/share/fonts/liberation/LiberationMono-Bold.ttf", 14)
    cw.add_grid(img, area_def, (10.0, 10.0), (2.0, 2.0), fontClr2, fill = "white", outline = 'lightblue',
                minor_outline = 'lightblue')

    # Save the image
    img.save(imgStr)

    # return the image object
    return img
# end definition

# definition-4
def imResize(imgStr):
    """
    What does this function do?
    Takes an image string and opens it in PIL, resizes the image and saves it with the same
    name and returns the image name as string again

    :param imgStr:
    :return: imgStr
    """
    import os, sys
    from PIL import Image

    img3 = Image.open(imgStr)
    baseWidth = 1024.0
    wPrcnt = (baseWidth / float(img3.size[0]))
    hSize = int((float(img3.size[1]) * float(wPrcnt)))
    webImg = img3.resize((int(baseWidth), hSize), Image.ANTIALIAS)
    webImg.save(imgStr)

    return imgStr
# end definition

# definition-5
def nc_write_sat_level_1_5(imgScn, outImgStr, prodStr):
    """
    What does this code do?
    This code takes the Satellite level 1.5 data as satpy-scene object from
    which the dask-array data is read and writes it to a netCDF file.

    :param imgScn:
    :param outImgStr:
    :param prodStr:
    :return:
    """
    # import necessary modules
    import dask
    import netCDF4 as nc4
    from datetime import datetime

    # Access the data from the scene object
    dataArray = imgScn[prodStr]
    satData = dataArray.compute()
    sData = satData.data
    lons_array, lats_array = satData.attrs['area'].get_lonlats()
    lon = lons_array[0]
    lat = lats_array[:,0]

    # Create the netCDF dataset template
    f = nc4.Dataset(outImgStr, 'w', format = 'NETCDF4')
    msg8_grp = f.createGroup('satellite_data')
    msg8_grp.description = "Satellite Level 1.5 Data for Channel: " + prodStr

    # Specify the dimensions
    # msg8_grp.createDimension('time', None)
    msg8_grp.createDimension('lons', len(lon))
    msg8_grp.createDimension('lats', len(lat))

    # Building the variables
    longitude = msg8_grp.createVariable('Longitude', 'f4', 'lons')
    latitude = msg8_grp.createVariable('Latitude', 'f4', 'lats')
    # startTime = msg8_grp.createVariables('Start_Time', 'i4', 'time')
    sat_data = msg8_grp.createVariable('Satellite_Channel_Data', 'f4', ('lats', 'lons'))
    # sat_data = msg8_grp.createVariable('Satellite_Channel_Data', 'f4', ('time', 'lons', 'lats'))

    # Ingest the actual data
    # longitude[:] = lons_array
    # latitude[:] = lats_array
    longitude[:] = lon
    latitude[:] = lat
    sat_data[:, :] = sData
    # sat_data[0,:,:] = sData
    # startTime = datetime.toordinal(satData.attrs['start_time'])

    # Adding Attributes
    longitude.units = 'degrees east'
    latitude.units = 'degrees north'
    # startTime.units = 'days since Jan 01, 0001'
    # From the Satellite Data
    sat_data.sensor = satData.attrs['sensor']
    sat_data.platform = satData.attrs['platform_name']
    sat_data.satellite_longitude = satData.attrs['satellite_longitude']
    sat_data.satellite_latitude = satData.attrs['satellite_latitude']
    sat_data.standard_name = satData.attrs['standard_name']
    sat_data.channel_name = satData.attrs['name']
    sat_data.units = satData.attrs['units']
    sat_data.wavelength = satData.attrs['wavelength']
    sat_data.resolution = satData.attrs['resolution']
    sat_data.fillValue = satData.attrs['_FillValue']
    sat_data.start_time = datetime.toordinal(satData.attrs['start_time'])
    sat_data.end_time = datetime.toordinal(satData.attrs['end_time'])
    sat_data.projection = satData.attrs['area'].proj_dict['proj']
    sat_data.ellipsoid = satData.attrs['area'].proj_dict['ellps']
    sat_data.upper_left_lon_lat = str(satData.attrs['area'].outer_boundary_corners[0])
    sat_data.upper_right_lon_lat = str(satData.attrs['area'].outer_boundary_corners[1])
    sat_data.lower_right_lon_lat = str(satData.attrs['area'].outer_boundary_corners[2])
    sat_data.lower_left_lon_lat = str(satData.attrs['area'].outer_boundary_corners[3])
    # Custom attributes
    sat_data.software = 'Pytroll_Satpy'
    sat_data.processed_by = 'MNRS, NCMRWF, Ministry of Earth Sciences'

    f.close()
# end definition

# definition-6
def deCmprss(msg1Src, exeDir, datDir, logDir, dateSnap, timeSnap):
    """
    What does this python function do?
    This function script decompresses the input HRIT files using EUMETSATś Public Wavelet Transform
    xRITDecompress and copies the uncompressed files to a temporary location for processing.
    ---MNRS---

    HISTORY:
    12-Feb-2018: Operational code defined as def()
    04-Jan-2018: Base python code for decompresion

    References:
    [1]. man xRITDecompress

    :param msg1Src: folder path of MSG-1 data for the said date & time
    :param exeDir: Folder path to xRIT decompression executable file
    :param datDir: folder path where the decompressed files are staged
    :param logDir: folder path where logs are present
    :param dateSnap: a date string in YYYYMMDD format; ex: 20180219
    :param timeSnap: a time string in hhmm format; ex: 2345
    :return:
    """
    #- Start
    # import necessary modules
    import os, sys, glob
    from datetime import datetime
    ## End of importing business

    print("\n \n \t \t \t STARTING THE DECOMPRESSION @ time: %s \t \t \t \n \n" % str(datetime.now()))

    ## Set directory paths
    srcDir = msg1Src

    deCmprssList = logDir + 'decompressedTimeSlots_' + dateSnap + '.txt'
    if not os.path.exists(deCmprssList):
        open(deCmprssList, 'w').close()
    # end if-not condition

    ## Go to exe directory
    os.chdir(exeDir)

    # Check if the time slot is already decompressed and available in datDir?
    if timeSnap in open(deCmprssList).read():
        print("\n...Skipping the decompression process for the time-stamp: %s" % timeSnap)
    else:
        print("\n...Decompressing for the time-stamp: %s" % timeSnap)
        print("\n \t \t Testing 123: \n \n ")
        searchComStr = srcDir, 'H-000-MSG1*' + timeSnap + '-C_'
        print(searchComStr)
        ## Start the loop business
        flist = glob.glob(os.path.join(srcDir, 'H-000-MSG1*' + timeSnap + '-C_'))

        for fname in flist:
            print('>>>Decompressing filename: %s' % fname)
            os.system('./xRITDecompress -s:%s' % fname)
            print('>>>Finished Decompresing filename: %s <<<' % fname)
        # end for-loop

        ## Move the decompressed files from exeDir to datDir
        mvCmd = "mv " + exeDir + "H-000*__*" + dateSnap + timeSnap + "-* " + datDir
        # mvCmd = "mv " + exeDir + "H-000*__*" + dateSnap +  "* " + datDir
        os.system(mvCmd)
        cpCmdEPI = "cp " + srcDir + "*EPI*" + dateSnap + timeSnap + "* " + datDir
        cpCmdPRO = "cp " + srcDir + "*PRO*" + dateSnap + timeSnap + "* " + datDir
        os.system(cpCmdEPI)
        os.system(cpCmdPRO)

        ## Mark the time stamp of decompressed files
        finDeTmsFile = logDir + "decompressedTimeSlots_" + dateSnap + ".txt"
        ff = open(finDeTmsFile, 'a+')
        ff.write("%s \n" % timeSnap)
        ff.close()
    # end if-condition to check whether the data is decompressed or not!
    print("deCmprss.py says: Finished with decompression of time-slot - %s - at Time: %s" % (timeSnap, str(datetime.now())))
# end-definition:q

# definition-7
def msg1Proc1_5(dateSnap, avail_times, fldrs):
    """
    What does this definition do?
    This script processes the raw MSG-1 Level 1.5 data to produces radiance/reflectance image
    files in netCDF,-4 geoTIFF & png file formats.

    :param dateSnap:
    :param avail_times:  A single string NOT an array
    :param fldrs:
    :return:
    """
    #- Start coding
    # import necessary modules
    import os, sys, glob
    #from satpy.utils import debug_on
    from satpy.scene import Scene
    from datetime import datetime
    from myDefinitions import nc_write_sat_level_1_5, embellish, imResize

    # Start the logic
    #debug_on()
    print("\n \t \t \t STARTING THE msg1Proc1_5 run @ time: %s \t \t \t \n \n" % str(datetime.now()))
    print("\n.Processing Date set is: %s" % dateSnap)

    #  Test whether all data folders are appropriately set or not.
    basDir, datDir, outDir, logDir, webDir, geoTdir, msg1Src, exeDir, GSHHS_ROOT, tmpDir = fldrs
    print("\n.Base directory is set to: %s" % basDir)
    print("\n.Data directory is set to %s" % datDir)
    print("\n.NetCDF output directory is set to: %s" % outDir)
    print("\n.Log directory is set to: %s" % logDir)
    print("\n.Web directory is set to: %s" % webDir)
    print("\n.GeoTiff directory is set to: %s" % geoTdir)
    print("\n.msg1Src directory is set to: %s" % msg1Src)
    print("\n.exeDir directory is set to: %s" % exeDir)
    print("\n.GSHHS directory is set to: %s" % GSHHS_ROOT)
    print("\n.tmpDir directory is set to: %s" % tmpDir)

    avail_times = str(avail_times).split()
    for tt in avail_times:
        try:
            # Start for-loop-1
            print("..Started processing for time: %s" % tt)
            searchStr = datDir + 'H-000-MSG1*' + dateSnap + tt + '-*'
            # searchStr = msg1Src + 'H-000-MSG1*' + dateSnap + tt + '*'

            files = glob.glob(searchStr)

            # Start reading filename in satpy
            scn = Scene(filenames=files, reader='hrit_msg')

            available_comps = scn.available_composite_names()
            channels_inverted = [s for s in available_comps if "_inv" in s]

            # add the remaining 3 non-inverted channels & 3d channel
            #allChnls = channels_inverted  + ["IR_016", "VIS006", "VIS008", "ir108_3d"]
            allChnls = channels_inverted  + ["ir108_3d"]

            # Save the individual channels (except HRV) as separate gray-scale GeoTIFF files..
            for ii in allChnls:
                try:
                    str(ii).split()
                    print("Working on channel: %s" % ii)
                    scn.load(str(ii).split())
                    indImg = scn.resample('India_SC')

                    # # Save as netCDF data
                    # outImgStr1 = outDir + 'ind_MSG1-Band_' + ii + '_' + dateSnap + '_' + tt + '.nc'
                    # nc_write_sat_level_1_5(indImg, outImgStr1, ii)

                    # # Save as Full Resolution GeoTIFF files
                    # outImgStr2 = geoTdir + 'ind_' + ii + '_' + dateSnap + '_' + tt + '.tiff'
                    # indImg.save_dataset(ii, filename = outImgStr2, writer = 'geotiff')

                    # Add graphics
                    # img2 = embellish(basDir, GSHHS_ROOT, outImgStr2, ii, dateSnap, tt)
                    # img2.save(outImgStr2)

                    # Save the data as resized png files
                    outImgStr3 = tmpDir + 'ind_' + ii + '_' + dateSnap + '_' + tt + '.png'
                    outImgStr3w = webDir + 'ind_' + ii + '_' + dateSnap + '_' + tt + '.png'
                    indImg.save_dataset(ii, filename = outImgStr3, writer = "simple_image")
                    outImgStr3 = imResize(outImgStr3)

                    # Add graphics
                    img3 = embellish(basDir, GSHHS_ROOT, outImgStr3, ii, dateSnap, tt)
                    img3.save(outImgStr3)

                    #  move the tmp files to proper web area
                    mv2WebCmd = 'mv ' + outImgStr3 + ' ' + outImgStr3w
                    os.system(mv2WebCmd)

                    # unload the read channel data
                    scn.unload(str(ii).split())
                    print("Finished processing for channel: %s " % ii)
                except:
                    print("Something went wrong with this Channel: %s" % ii)
                    continue
                # end try-except block
            #end for-loop
            print("Finished processing for time-stamp: %s" % tt)
        except:
            print("Something went wrong with this time: %s" % tt)
            continue
        # end try-except block
    # end for-loop
# end-definition

# definition-8
def msg1RGBProc(dateSnap, avail_times, fldrs):
    """
    What does this definition do?
    This script processes the raw MSG-1 data into RGB Data Products in netCDF-4, geoTIFF &
    png file formats

    :param dateSnap:
    :param avail_times: A single string NOT an array
    :param fldrs:
    :return:
    """
    #-Start coding
    # start the logic
    import os, sys, glob
    #from satpy.utils import debug_on
    from satpy.scene import Scene
    from datetime import datetime
    from myDefinitions import nc_write_sat_level_1_5, embellish, imResize

    # Start the logic
    #debug_on()
    print("\n \t \t \t STARTING THE msg1RGBProc run @ time: %s \t \t \t \n \n" % str(datetime.now()))
    print("\n.Processing Date set is: %s" % dateSnap)

    #  Test whether all data folders are appropriately set or not.
    basDir, datDir, outDir, logDir, webDir, geoTdir, msg1Src, exeDir, GSHHS_ROOT, tmpDir = fldrs
    print("\n.Base directory is set to: %s" % basDir)
    print("\n.Data directory is set to %s" % datDir)
    print("\n.NetCDF output directory is set to: %s" % outDir)
    print("\n.Log directory is set to: %s" % logDir)
    print("\n.Web directory is set to: %s" % webDir)
    print("\n.GeoTiff directory is set to: %s" % geoTdir)
    print("\n.msg1Src directory is set to: %s" % msg1Src)
    print("\n.exeDir directory is set to: %s" % exeDir)
    print("\n.GSHHS directory is set to: %s" % GSHHS_ROOT)
    print("\n.tmpDir directory is set to: %s" % tmpDir)

    avail_times = str(avail_times).split()
    for tt in avail_times:
        # Start for-loop-1
        print("..Started processing for time: %s" % tt)
        searchStr = datDir + 'H-000-MSG1*' + dateSnap + tt + '-*'
        files = glob.glob(searchStr)

        # Start reading filename in satpy
        scn = Scene(filenames=files, reader='hrit_msg')

        #for composite in ['airmass', 'ash', 'cloudtop', 'colorized_ir_clouds', 'convection', 'day_microphysics', 'dust', 'fog', 'night_fog', 'natural_color', 'ir_overview']:
        for composite in ['fog','ash', 'cloudtop', 'colorized_ir_clouds', 'day_microphysics']:
            if (composite == 'airmass'):
                prodStr = 'airM'
                capStr = 'Air Mass'
            elif (composite == 'ash'):
                prodStr = 'ASH'
                capStr = 'ASH'
            elif (composite == 'cloudtop'):
                prodStr = 'CTOP'
                capStr = 'Cloud Top View'
            elif (composite == 'cloud_optical_thickness'):
                prodStr = 'COP'  # 2 much Problematic
                capStr = 'Cloud Optical Thickness'
            elif (composite == 'cloud_top_temperature'):
                prodStr = 'CTT'  # problematic
                capStr = 'Cloud Top Temperature'
            elif (composite == 'cloud_top_height'):
                prodStr = 'CTH'
                capStr = 'Cloud Top Height'
            elif (composite == 'cloudtype'):
                prodStr = 'CType'
                capStr = 'Cloud Type'
            elif (composite == 'cloud_top_pressure'):
                prodStr = 'CTP'
                capStr = 'Cloud Top Pressure'
            elif (composite == 'cloud_top_phase'):
                prodStr = 'CTPh'
                capStr = 'Cloud Top Phase'
            elif (composite == 'cloudmask'):
                prodStr = 'CMask'
                capStr = 'Cloud Mask'
            elif (composite == 'colorized_ir_clouds'):
                prodStr = 'C-IR'
                capStr = 'Colourised IR Clouds'
            elif (composite == 'convection'):
                prodStr = 'CON'
                capStr = 'Convection Activity'
            elif (composite == 'day_microphysics'):
                prodStr = 'd-MicPhy'
                capStr = 'Day Microphysics'
            elif (composite == 'dust'):
                prodStr = 'dust'
                capStr = 'DUST'
            elif (composite == 'fog'):
                prodStr = 'FOG'
                capStr = 'Fog Activity'
            elif (composite == 'night_fog'):
                prodStr = 'NFog'
                capStr = 'Night Time Fog Activity'
            elif (composite == 'natural_color'):
                prodStr = 'NAT'
                capStr = 'Quasi True Colour'
            elif (composite == 'ir_overview'):
                prodStr = 'IR'
                capStr = 'False Colour Composite of IR'
            # end if condition

            try:
                # Load the scene
                scn.load([composite])

                # India Specific Scene
                indScn = scn.resample("India_SC")
                indScn.load([composite])

                # # Save as netCDF data ---- TO BE IMPLEMENTED ----
                # outImgStr1 = outDir + 'ind_MSG-1_RGB_' + prodStr + '_' + dateSnap + '_' + tt + '.nc'
                # # indImg.save_datasets(writer = 'cf', filename = outImgStr1)
                # nc_write_sat_level_1_5(indScn, outImgStr1, prodStr)

                # # Save as Full Resolution GeoTIFF files
                # outImgStr2 = geoTdir + 'ind_MSG-1_RGB_' + prodStr + '_' + dateSnap + '_' + tt + '.tiff'
                # indScn.save_dataset(composite, filename = outImgStr2, writer = 'geotiff')
                # # Add graphics
                # # img2 = embellish(basDir, GSHHS_ROOT, outImgStr2, capStr, dateSnap, tt)
                # # img2.save(outImgStr2)

                # Save the data as resized png files
                outImgStr3 = tmpDir + 'ind_MSG1_RGB_' + prodStr + '_' + dateSnap + '_' + tt + '.png'
                outImgStr3w= webDir + 'ind_MSG1_RGB_' + prodStr + '_' + dateSnap + '_' + tt + '.png'
                indScn.save_dataset(composite, filename = outImgStr3, writer = "simple_image")
                outImgStr3 = imResize(outImgStr3)
                # Add graphics
                img3 = embellish(basDir, GSHHS_ROOT, outImgStr3, capStr, dateSnap, tt)
                img3.save(outImgStr3)

                #  move the tmp files to proper web area
                mv2WebCmd = 'mv ' + outImgStr3 + ' ' + outImgStr3w
                os.system(mv2WebCmd)

                # unload the read channel data
                scn.unload([composite])
                indScn.unload([composite])
                print("Finished processing for RGB Composite: %s " % composite)
            except:
                print("Something went wrong with this RGB composite: %s" % composite)
                continue
            # end try-except block
        # end for-loop

        #Finished time slots
        finTmStmps = [tt]
        print("\n.Reading Finished Time slots as: %s" % finTmStmps)
        finTmsFile = logDir + "finishedTimeSlots_" + dateSnap + ".txt"
        fp = open(finTmsFile, 'a+')
        for item in finTmStmps:
            fp.write("%s \n" % item)
        # end for loop to write
        fp.close()
        # end for-loop
    # end-for-loop
    print("msg1RGBProc() says: Finished with processing of time-slot - %s - at: %s " % (tt, str(datetime.now())))
# end-definition

# definition-9
def msg1NDVI(dateSnap, avail_times, fldrs):
    """
    What does this function do?
    This definition/function is meant for computing NDVI from SEVIRI data

    Ref: https://nbviewer.jupyter.org/github/pytroll/pytroll-examples/blob/master/satpy/hrit_msg_tutorial.ipynb

    :param dateSnap:
    :param avail_times:
    :param fldrs:
    :return: NDVI
    """

    # Start the logic
    import os, sys, glob
    #from satpy.utils import debug_on
    from satpy.scene import Scene
    from satpy.dataset import combine_metadata
    from datetime import datetime
    from trollimage.colormap import greys, greens
    from trollimage.image import Image
    from myDefinitions import nc_write_sat_level_2, embellish, imResize

    #debug_on()

    print("\n \t \t \t STARTING THE msg1NDVI run @ time: %s \t \t \t \n \n" % str(datetime.now()))
    print("\n.Processing Date set is: %s" % dateSnap)

    #  Test whether all data folders are appropriately set or not.
    basDir, datDir, outDir, logDir, webDir, geoTdir, msg1Src, exeDir, GSHHS_ROOT, tmpDir = fldrs
    print("\n.Base directory is set to: %s" % basDir)
    print("\n.Data directory is set to %s" % datDir)
    print("\n.NetCDF output directory is set to: %s" % outDir)
    print("\n.Log directory is set to: %s" % logDir)
    print("\n.Web directory is set to: %s" % webDir)
    print("\n.GeoTiff directory is set to: %s" % geoTdir)
    print("\n.msg1Src directory is set to: %s" % msg1Src)
    print("\n.exeDir directory is set to: %s" % exeDir)
    print("\n.GSHHS directory is set to: %s" % GSHHS_ROOT)
    print("\n.tmpDir directory is set to: %s" % tmpDir)

    avail_times = str(avail_times).split()
    for tt in avail_times:
        # Start for-loop-1
        print("..Started processing for time: %s" % tt)
        searchStr = datDir + 'H-000-MSG1*' + dateSnap + tt + '-*'
        print("\n \t \t Testing 123: \n \n ")
        print(searchStr)
        files = glob.glob(searchStr)
        #print("\n Testing 123: \n")
        #print(files)

        # Start reading filename in satpy
        scn = Scene(filenames=files, reader='hrit_msg')

        #  start the NDVI computation
        scn.load(['VIS006', 0.6])
        scn.load(['VIS008', 0.8])
        ndvi = (scn[0.8] - scn[0.6]) / (scn[0.8] + scn[0.6])
        ndvi.attrs = combine_metadata(scn[0.6], scn[0.8])
        scn['ndvi'] = ndvi

        composite = 'ndvi'
        prodStr = 'NDVI'
        capStr = 'NDVI'

        # resample the data to Indian region
        indScn = scn.resample('India_SC')

        #  save the data
        # # # Save as netCDF data ---- TO BE IMPLEMENTED ----
        # outImgStr1 = outDir + 'ind_MSG-1_RGB_' + prodStr + '_' + dateSnap + '_' + tt + '.nc'
        # nc_write_sat_level_2(indScn, outImgStr1, prodStr)
        #
        # # Save as Full Resolution GeoTIFF files
        # outImgStr2 = geoTdir + 'ind_MSG-1_RGB_' + prodStr + '_' + dateSnap + '_' + tt + '.tiff'
        # indScn.save_dataset(composite, filename = outImgStr2, writer = 'geotiff')
        # # Add graphics
        # # img2 = embellish(basDir, GSHHS_ROOT, outImgStr2, capStr, dateSnap, tt)
        # # img2.save(outImgStr2)

        # Save the data as resized png files
        outImgStr3 = tmpDir + 'ind_MSG1_RGB_' + prodStr + '_' + dateSnap + '_' + tt + '.png'
        outImgStr3w = webDir + 'ind_MSG1_RGB_' + prodStr + '_' + dateSnap + '_' + tt + '.png'

        #  Apply color palette from trollimage
        ndvi_data = indScn['ndvi'].compute().data
        ndvi_img = Image(ndvi_data, mode = "L")
        # greys.set_range(ndvi_data.min(), -0.00001)
        # greens.set_range(0,ndvi_data.max())
        greys.set_range(-0.8, -0.00001)
        greens.set_range(0, 0.8)
        my_cm = greys + greens
        ndvi_img.colorize(my_cm)
        ndvi_img.save(outImgStr3)
        # indScn.save_dataset(composite, filename = outImgStr3, writer = "simple_image")
        outImgStr3 = imResize(outImgStr3)
        # Add graphics
        img3 = embellish(basDir, GSHHS_ROOT, outImgStr3, capStr, dateSnap, tt)
        img3.save(outImgStr3)

        #  move the tmp files to proper web area
        mv2WebCmd = 'mv ' + outImgStr3 + ' ' + outImgStr3w
        os.system(mv2WebCmd)

        print("msg1NDVI() says: Finished with processing of time-slot - %s - at: %s " % (tt, str(datetime.now())))
    # end for-loop
# end definition

#  definition-10
def nc_write_sat_level_2(imgScn, outImgStr, prodStr):
    """
    What does this code do?
    This code takes the Satellite level 2 derived data as satpy-scene object from
    which the dask-array data is read and writes it to a netCDF file.

    :param imgScn:
    :param outImgStr:
    :param prodStr:
    :return:
    """
    # import necessary modules
    import dask
    import netCDF4 as nc4
    from datetime import datetime

    # Access the data from the scene object
    dataArray = imgScn[prodStr]
    satData = dataArray.compute()
    sData = satData.data
    lons_array, lats_array = satData.attrs['area'].get_lonlats()
    lon = lons_array[0]
    lat = lats_array[:,0]

    # Create the netCDF dataset template
    f = nc4.Dataset(outImgStr, 'w', format = 'NETCDF4')
    msg8_grp = f.createGroup('satellite_data')
    msg8_grp.description = "Satellite Level 2 Derived Data Product: " + prodStr

    # Specify the dimensions
    # msg8_grp.createDimension('time', None)
    msg8_grp.createDimension('lons', len(lon))
    msg8_grp.createDimension('lats', len(lat))

    # Building the variables
    if prodStr == 'ndvi':
        longitude = msg8_grp.createVariable('Longitude', 'f4', 'lons')
        latitude = msg8_grp.createVariable('Latitude', 'f4', 'lats')
        # startTime = msg8_grp.createVariables('Start_Time', 'i4', 'time')
        sat_data = msg8_grp.createVariable('Satellite_Channel_Data', 'f4', ('lats', 'lons'))
        # sat_data = msg8_grp.createVariable('Satellite_Channel_Data', 'f4', ('time', 'lons', 'lats'))

        # Ingest the actual data
        # longitude[:] = lons_array
        # latitude[:] = lats_array
        longitude[:] = lon
        latitude[:] = lat
        sat_data[:, :] = sData
        # sat_data[0,:,:] = sData
        # startTime = datetime.toordinal(satData.attrs['start_time'])

        # Adding Attributes
        longitude.units = 'degrees east'
        latitude.units = 'degrees north'
        # startTime.units = 'days since Jan 01, 0001'
        # From the Satellite Data
        sat_data.sensor = satData.attrs['sensor']
        sat_data.platform = satData.attrs['platform_name']
        sat_data.satellite_longitude = satData.attrs['satellite_longitude']
        sat_data.satellite_latitude = satData.attrs['satellite_latitude']
        sat_data.standard_name = 'normalized_difference_vegetation_index'
        sat_data.product_name = satData.attrs['name']
        sat_data.units = 'ratio'
        sat_data.resolution = satData.attrs['resolution']
        sat_data.fillValue = satData.attrs['_FillValue']
        sat_data.start_time = datetime.toordinal(satData.attrs['start_time'])
        sat_data.end_time = datetime.toordinal(satData.attrs['end_time'])
        sat_data.projection = satData.attrs['area'].proj_dict['proj']
        sat_data.ellipsoid = satData.attrs['area'].proj_dict['ellps']
        sat_data.upper_left_lon_lat = str(satData.attrs['area'].outer_boundary_corners[0])
        sat_data.upper_right_lon_lat = str(satData.attrs['area'].outer_boundary_corners[1])
        sat_data.lower_right_lon_lat = str(satData.attrs['area'].outer_boundary_corners[2])
        sat_data.lower_left_lon_lat = str(satData.attrs['area'].outer_boundary_corners[3])
        # Custom attributes
        sat_data.software = 'Pytroll_Satpy'
        sat_data.processed_by = 'MNRS, NCMRWF, Ministry of Earth Sciences'
# end if-condition
    f.close()
    print("\n Finished writing the netCDF file for NDVI product")
# end definition

# EOF
