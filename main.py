#!/usr/bin/env python
__author__ = "raghav"

'''
What does this code do?
This program is the main calling program for processing MSG data

HISTORY:
09-Jul-2018: Migration to Mihir
14-Mar-2018: subprocess based version
09-Feb-2018: Operational code for main
07-Jan-2018: Base code for msg1Read.py

References:
[1]. Logging: https://stackoverflow.com/questions/11124093/redirect-python-print-output-to-logger?rq=1

'''
#- Start code
# Start the main definition
def main():
    # import necessary system modules
    import os, glob
    import numpy as np
    from datetime import datetime
    # import necessary custom made def modules
    from myDefinitions import dateTimeManipulation, mkFldrs, deCmprss, msg1Proc1_5, msg1NDVI, msg1RGBProc
    # end of module importing business

    # print stuff for log header
    print("\n \n \t \t \t STARTING THE MAIN RUN @ time: %s \t \t \t \n \n" % str(datetime.now()))

    # manipulate datetime as per needs
    dateSnap, dt, yr, mm, dd, tm, hrs, mins, dt_2, yr_2, mm_2, dd_2, tm_2, hrs_2, mins_2 = dateTimeManipulation()

    # Set-up the folders & log files
    basDir, datDir, outDir, logDir, webDir, geoTdir, msg1Src, exeDir, GSHHS_ROOT = mkFldrs(dateSnap)
    fldrs = [basDir, datDir, outDir, logDir, webDir, geoTdir, GSHHS_ROOT]

    # Generate the time stamps when MSG1 data is available
    msgTmStmps = [str(datetime(int(yr),int(mm),int(dd),hr, mn).time().strftime("%H%M")) for hr in range(0,24) for mn in range(0,60,15)]
    msgTmStmps_int = np.array(map(int, msgTmStmps)) # Just maintain the same as integers for easier comparison
    finTmStmps = []

    # Get unprocessed data time-stamps in the last 2 hour
    pastTmStmp = int(hrs_2+mins_2)
    presentTmStmp = int(hrs+mins)
    indx = np.array(np.where((msgTmStmps_int >= pastTmStmp) & (msgTmStmps_int <= presentTmStmp))).squeeze()
    procTmStmps = msgTmStmps[indx[0]:indx[-1]]

    # Start processing for Level-1.5 radiance/reflectance data
    for tt in procTmStmps:
        print("\n > Working on following time-stamps: %s " % procTmStmps)
        print("\n >> Searching for time-stamp: *%s*" % tt)
        # Check if the time slot is already processed?
        if tt in open(logDir + 'finishedTimeSlots_' + dateSnap + '.txt').read():
            print("\n >>> Already processed.. Skipping the time-step: %s" % tt)
        else:
            # Check if both EPI and PRO are available
            print("\n >>>> Processing for time: %s" % tt)
            # Check if PRO files exist for the time-slot or not
            if (len([f for f in glob.iglob(msg1Src + '*PRO*' + tt + '-*')])):
                print("\n >>>>> PRO file exists, Lemme check for EPI now")
               # Check if EPI files exists for the time-slot or not
                if (len([f for f in glob.iglob(msg1Src + '*EPI*' + tt + '-*')])):
                    print("\n >>>>> EPI file too exists, Lemme proceed with decompression & processing \n \n")
                    # decompress the files first
                    deCmprss(msg1Src, exeDir, datDir, logDir, dateSnap, tt)
                    os.chdir(basDir)
                    print("\n <<<<< Finished decompressing scheme, exiting from deCmprss() \n")
                    # call the main processing modules.
                    msg1Proc1_5(dateSnap, tt, fldrs)
                    msg1NDVI(dateSnap, tt, fldrs)
                    msg1RGBProc(dateSnap,tt, fldrs)
                    print("\n <<<<< Finished processing for RGB Data products, exiting from msg1RGBProc() \n")
                else:
                    print("\n >>>> EPI File does NOT exist, Skipping the processing for this time step: %s" % tt)
                # end if-EPI
            else:
                print("\n >>>> PRO File does NOT exist, Skipping the EPI Search & Processing altogether for time step: %s" % tt)
            # end if-PRO
        # end if-condition
        print("\n < Finished working with time-stamp: %s" % tt)
        os.chdir(basDir)
    # end-for loop
#  end-main

if __name__ == "__main__":
    from datetime import datetime
    main()
    print("\n \n \t \t \t FINISHED THE MAIN RUN @ time: %s \t \t \t \n \n" % str(datetime.now()))

###EOF###

