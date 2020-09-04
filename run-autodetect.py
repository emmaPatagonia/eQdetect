import os, glob
import numpy as np
from copy import deepcopy as cp
from obspy.core import read, UTCDateTime, Stream
from scripts import autodetect, runHypo71
#from pygema.read import get_stations_info, get_waveforms

class bcolors:
  HEADER = '\033[95m'
  OKBLUE = '\033[94m'
  OKGREEN = '\033[92m'
  WARNING = '\033[93m'
  FAIL = '\033[91m'
  ENDC = '\033[0m'
  BOLD = '\033[1m'
  UNDERLINE = '\033[4m'


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# READ WAVEFORMS
#starttime = UTCDateTime("2019-09-29 11:30:00")
#endtime   = UTCDateTime("2019-10-01 11:30:00")
starttime = UTCDateTime("2019-09-30 18:45:00") 
endtime   = UTCDateTime("2019-09-30 19:25:00")

msfile = "src/msfiles/STREAM_2019.09.29_2019.10.01"


st = read(msfile, starttime=starttime, endtime=endtime )


# READ STATION INFORMATION
network_info = np.loadtxt("src/stations.net", dtype="str")
networks = network_info[0]
stations = network_info[1]
#st = Stream() 
#for network,station in zip(networks,stations):
#  print("    -> %s.%s    %s   %s" % (network,station,starttime.strftime("%Y-%m-%d %H:%M:%S"),endtime.strftime("%Y-%m-%d %H:%M:%S")) )
#  st += get_waveforms(network,station,starttime,endtime)


# SET PARAMETERS FOR PRE-PROCESSING OF RAWDATA
freqmin = 3.5
freqmax = 10
tapering = 0.05



# SET PARAMETERS FOR TRIGGER COINCIDENT
time_window_length = 30*60

sta = 0.5
lta = 10
thr_on = 3.
thr_off = 1.5

min_num_stations = 4
deadtime_between_coincidences = 10


# SET START/END TIMES WHEN EVENT IS FOUND AND LOCALIZE USING HYPO71
time_before = cp(lta) # time before the trigger coincidence when trimming the short segment
time_after = 60*2     # time after the trigger coincidence when trimming the short segment
deadtime_after_pphase = 3 # deadtime after the trigger coincidence when trimming for s-phase


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  1) RUN TRIGGER COINCIDENT IN A LONG TIME-WINDOW (take in consideration the lta parameter length)
#  2) IF ANY COINCIDENCE EXISTS, COMPUTE STA/LTA FOR P-PHASE AND S-PHASE IN A SHORT TIME WINDOW
#  3) THEN, RETURN DICTIONARY OF EVENTS IN HYPO71 FORMAT

coincidences_dict = autodetect(st, freqmin, freqmax, tapering, sta, lta, thr_on, thr_off, min_num_stations, deadtime_between_coincidences, time_before, time_after, deadtime_after_pphase)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# LOCALIZAR EVENTOS DETECTADOS USANDO HYPO71 

# RUN HYPO71
if len(coincidences_dict)>0:
  events_dict = runHypo71(coincidences_dict, maxgap=360)

  # LOOP OVER EACH EVENT
  for evid in events_dict:
    infodict = events_dict[evid]

    evtime = infodict['origin_time']
    evlon  = infodict['evlon']
    evlat  = infodict['evlat']
    evdep  = infodict['evdep']
    evmag  = -99 
    evnstats = infodict['nphases']
    evgap  = infodict['gap']
    evrms  = infodict['rms']
    everrx = -99
    everry = -99
    everrz = -99
    status = "automatic"
    #print(bcolors.WARNING+"[insert]  %s  %.4f %.4f  %.1f km   Ml %.1f  %i %.1f %.1f    %.1f km %.1f km %.1f km  (%s) " % (evtime.strftime("%Y-%m-%d %H:%M:%S"), evlon, evlat, evdep, evmag, evnstats, evgap, evrms, everrx, everry, everrz, status  ) + bcolors.ENDC )



