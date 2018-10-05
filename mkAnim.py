#!/usr/bin/env python

"""
What does this script do?


"""

# import necessary custom made def modules
import numpy as np
import imageio
from datetime import datetime
from myDefinitions import dateTimeManipulation, mkFldrs
import os, glob, sys

# Start the logic
dateSnap, dt, yr, mm, dd, tm, hrs, mins, dt_2, yr_2, mm_2, dd_2, tm_2, hrs_2, mins_2 = dateTimeManipulation()
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

with imageio.get_writer('/home/sateeshm/myPY/animation/MSG_animation.gif', mode='I', fps=5) as writer:
    for counter in range(len(procTmStmps)):
        for filename in glob.glob(dir_path + '*RGB_IR*' + dateSnap + '_' + procTmStmps[counter] + '.png'):
            image = imageio.imread(filename)
            writer.append_data(image)
        # end for-loop
    # end for-loop
# end with-loop

# EOF
