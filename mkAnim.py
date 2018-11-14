#!/usr/bin/env python
__author__ = "sateeshm", "raghav"

"""
What does this script do?
This script was written to create GIF animations from the live data. The expected
output would be a simple GIF file and/or a mp4 movie data.
This script is intended to be the base for web-query based on-demand visualization in 
near future.

HISTORY:
04-Oct-2018: First prototype version initiated by M. Sateesh.
05-Oct-2018: Upgradation and improvements by MNRS for possible generalization to 
on-the-fly web-based requests.

References:
[1].  https://stackoverflow.com/questions/753190/programmatically-generate-video-or-animated-gif-in-python
[2]. 

"""

# import necessary custom made def modules
import numpy as np
import imageio
from datetime import datetime
import os, glob, sys

# Start the logic
dateSnap, dt, yr, mm, dd, tm, hrs, mins, dt_2, yr_2, mm_2, dd_2, tm_2, hrs_2, mins_2 = set_date_time_stamps()
dir_path = '/home/raghav/Work/satViz/web/' + dateSnap + '/'

# Generate the time stamps when MSG1 data is available
msgTmStmps = [str(datetime(int(yr), int(mm), int(dd), hr, mn).time().strftime("%H%M")) for hr in range(0, 24) for mn in
              range(0, 60, 15)]
msgTmStmps_int = np.array(map(int, msgTmStmps))  # Just maintain the same as integers for easier comparison

# Get unprocessed data time-stamps in the last 2 hour
pastTmStmp = int(hrs_2 + mins_2)
presentTmStmp = int(hrs + mins)
indx = np.array(np.where((msgTmStmps_int >= pastTmStmp) & (msgTmStmps_int <= presentTmStmp))).squeeze()
procTmStmps = msgTmStmps[indx[0]:indx[-1]]

# build the gif animation

with imageio.get_writer('/home/raghav/Work/satViz/animations/MSG_animation.gif', mode='I', fps=4) as writer:
    for counter in range(len(procTmStmps)):
        for filename in glob.glob(dir_path + '*RGB_CON**' + dateSnap + '_' + procTmStmps[counter] + '.png'):
            image = imageio.imread(filename)
            writer.append_data(image)
        # end for-loop
    # end for-loop
# end with-loop

# Local definitions
def set_date_time_stamps(start_time, end_time, prodStr):
    """
    What does this definition do?
    Manipulates date time strings as per needs

    :param: start_time: a date time string in YYYYMMDDHHmm format, example: 201810011345
    :param: end_time: a date time string same as start_time format, example: 201810011830
    :param: prodStr: The product that you are trying to create the animation for, example: CON
            it should conform to the known prodStr.
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

    # Past 2 hour time
    dt_2 = str(datetime.utcnow() - timedelta(hours = 1.0)).split()
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


# EOF
