import os, glob
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy as cp
from obspy.core import read, UTCDateTime, Stream
from obspy.signal import trigger 
from obspy.geodetics.base import degrees2kilometers, calc_vincenty_inverse
from hypo71 import hypo71

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

class bcolors:
  HEADER = '\033[95m'
  OKBLUE = '\033[94m'
  OKGREEN = '\033[92m'
  WARNING = '\033[93m'
  FAIL = '\033[91m'
  ENDC = '\033[0m'
  BOLD = '\033[1m'
  UNDERLINE = '\033[4m'


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def autodetect(st, freqmin, freqmax, tapering, sta, lta, thr_on, thr_off, min_num_stations, deadtime_between_coincidences, time_before, time_after, deadtime_after_pphase):
  """ 
  + RUN TRIGGER COINCIDENT IN A LONG TIME-WINDOW (take in consideration the lta parameter length)
  + IF ANY COINCIDENCE EXISTS, COMPUTE STA/LTA FOR P-PHASE AND S-PHASE IN A SHORT TIME WINDOW
  + THEN, RETURN DICTIONARY OF EVENTS IN HYPO71 FORMAT
  """

  # PRE-PROCESSING OF RAWDATA
  print(bcolors.BOLD + "\n+ Pre-processing rawdata..." + bcolors.ENDC)
  st.detrend("demean")
  st.taper(max_percentage=tapering, type="hann")
  st.merge(method=1, fill_value='interpolate')
  st.filter("bandpass",freqmin=freqmin,freqmax=freqmax)
  st.sort()


  # REMOVE TAPERED CORNERS
  for tr in st:
    dt = tapering*(tr.stats.endtime-tr.stats.starttime)
    tr.trim(tr.stats.starttime+dt, tr.stats.endtime-dt)


  # RUN COINCIDENCE TRIGGER (STA/LTA + CORRELATION)
  print(bcolors.BOLD + "\n+ Running trigger coincidence..." + bcolors.ENDC)
  st_z = st.copy().select(channel="*Z")
  output = trigger.coincidence_trigger("recstalta", thr_on=thr_on, thr_off=thr_off, 
                                        stream=st_z, 
                                        thr_coincidence_sum=min_num_stations, sta=sta, lta=lta,  
                                        trigger_off_extension=0, similarity_threshold=0.7, 
                                        details=True)

  # CREATE OUTPUT DICTIONARY OF EVENTS
  coincidences_dict = {}

  # LOOP OVER EACH EVENT (if anyone is found)
  timeX_old = 0
  evnum = 1
  for event in output:
    timeX = event['time']
    traces = event['stations']


    # IF EVENT IS REPEATED, THEN CONTINUE (JUMP) TO THE NEXT ITERATION
    if abs(timeX - timeX_old) < deadtime_between_coincidences:
      continue
    timeX_old = cp(timeX)

    # PRINT TRIGGER COINCIDENCE TIME
    print( bcolors.BOLD + "\n[%i] Coincident found at %s ..." % (evnum, timeX.strftime("%Y-%m-%d %H:%M:%S")) + bcolors.ENDC   )


    # LOOP OVER EACH STATION + CUT SHORT SEGMENT
    picks_list = []
    for stname in traces:
      try:
        st_short = st.copy().select(station=stname)
        t1 = timeX - time_before
        t2 = timeX + time_after
        st_short.trim(t1, t2)

        ######## RUN STA/LTA FOR P-PHASE ########
        tr_z = st_short.select(channel="*Z")[0]
        cft = trigger.recursive_sta_lta(tr_z.data, int(tr_z.stats.sampling_rate*sta), int(tr_z.stats.sampling_rate*lta))
        on_off = trigger.trigger_onset(cft, thr_on, thr_off)
        #trigger.plot_trigger(tr_z, cft,thr_on, thr_off, show=True)
        if len(on_off) > 0:
          weight = 0
          if len(on_off) > 1: 
            weight = 1
            #print(bcolors.WARNING + "(warning: more than one alert for P-phase at %s)"%(tr_z.stats.station) + bcolors.ENDC  )

          # DEFINE ARRIVAL TIME OF THE P-PHASE AS THE FIRST TRIGGER FOUND
          alert_P = tr_z.times("UTCDateTime")[on_off[0][0]]   
          alert_pattern = "%s P %.2f %i" % (tr_z.stats.station, alert_P.timestamp, weight)



          ######## RUN STA/LTA FOR S-PHASE (north component) ######## 
          flag_Snorth = False
          try:
            tr_n = st_short.select(channel="*N")[0]
            tr_n.trim(alert_P+deadtime_after_pphase, tr_n.stats.endtime)
            cft_n = trigger.recursive_sta_lta(tr_n.data, int(tr_n.stats.sampling_rate*sta), int(tr_n.stats.sampling_rate*lta))
            on_off_Sn = trigger.trigger_onset(cft_n, thr_on, thr_off)
            #trigger.plot_trigger(tr_n, cft_n,thr_on, thr_off, show=True)
            if len(on_off_Sn) > 0:
              flag_Snorth = True
              weight_Sn = 0
              if len(on_off_Sn) > 1: 
                weight_Sn = 1
                #print(bcolors.WARNING + "(Warning: more than one alert for S-phase north component)" + bcolors.ENDC)

              # DEFINE ARRIVAL TIME OF THE S-PHASE AS THE FIRST TRIGGER FOUND
              alert_Sn = tr_n.times("UTCDateTime")[on_off_Sn[0][0]]   

          except:
            continue


          ######## RUN STA/LTA FOR S-PHASE (east component) ########
          flag_Seast = False
          try:
            tr_e = st_short.select(channel="*E")[0]
            tr_e.trim(alert_P+deadtime_after_pphase, tr_e.stats.endtime)
            cft_e = trigger.recursive_sta_lta(tr_e.data, int(tr_e.stats.sampling_rate*sta), int(tr_e.stats.sampling_rate*lta))
            on_off_Se = trigger.trigger_onset(cft_e, thr_on, thr_off)
            #trigger.plot_trigger(tr_e, cft_e,thr_on, thr_off, show=True)
            if len(on_off_Se) > 0:
              flag_Seast = True
              weight_Se = 0
              if len(on_off_Se) > 1: 
                weight_Se = 1
                #print(bcolors.WARNING + "(Warning: more than one alert for S-phase east component)" + bcolors.ENDC)

              # DEFINE ARRIVAL TIME OF THE S-PHASE AS THE FIRST TRIGGER FOUND
              alert_Se = tr_e.times("UTCDateTime")[on_off_Se[0][0]]   

          except:
            continue


          # DEFINE ARRIVAL TIME OF THE S-PHASE BETWEEN THE NORTH AND EAST COMPONENT
          if flag_Snorth and flag_Seast:
            if weight_Sn<weight_Se:
              if alert_Sn-alert_P<=deadtime_after_pphase and alert_Sn-alert_P>0:
                alert_S = cp(alert_Sn)
                weight_S = cp(weight_Sn)
                alert_pattern += " S %.2f %i" % (alert_S.timestamp, weight_S)      

            elif weight_Sn>weight_Se:
              if alert_Se-alert_P<=deadtime_after_pphase and alert_Se-alert_P>0:
                alert_S = cp(alert_Se)
                weight_S = cp(weight_Se)
                alert_pattern += " S %.2f %i" % (alert_S.timestamp, weight_S)      

            else:
              timediff = abs(alert_Sn-alert_Se)
              if timediff <= 3:
                if alert_Sn <= alert_Se:
                  alert_S = alert_Sn+timediff/2.
                else:
                  alert_S = alert_Sn-timediff/2.
                weight_S = cp(weight_Sn)
                alert_pattern += " S %.2f %i" % (alert_S.timestamp, weight_S)      

          elif flag_Snorth and not flag_Seast:
            if alert_Sn-alert_P<=deadtime_after_pphase and alert_Sn-alert_P>0:
              alert_S = cp(alert_Sn)
              alert_pattern += " S %.2f %i" % (alert_S.timestamp, weight_Sn)      

          elif not flag_Snorth and flag_Seast:
            if alert_Se-alert_P<=deadtime_after_pphase and alert_Se-alert_P>0:
              alert_S = cp(alert_Se)
              alert_pattern += " S %.2f %i" % (alert_S.timestamp, weight_Se)      


          # PRINT OUTPUT
          print(" "*4 + bcolors.OKGREEN + alert_pattern + bcolors.ENDC)
          picks_list.append(alert_pattern)


      except:
        continue

    # APPEND PICKS TO DICTIONARY
    evname = "Event_%03i" % (evnum)
    coincidences_dict.update({evname: cp(picks_list)})
    evnum += 1

  return coincidences_dict




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


def export_picksfile(coincidences_dict, pickfile="picks.txt"):

  if len(coincidences_dict)>0:

    if os.path.isfile(pickfile):
      os.remove(pickfile)

    outfile = open(pickfile, "w")
    for event in coincidences_dict:
      outfile.write("#"*40)
      outfile.write("\n")

      for alert_pattern in coincidences_dict[event]:
        alert_pattern += "\n"
        outfile.write(alert_pattern)

    outfile.write("\n")
    outfile.close()

    print(bcolors.BOLD + "\n%+ i events found ... exporting %s file\n" % (len(coincidences_dict), pickfile) + bcolors.ENDC)

  else:
    print(bcolors.BOLD + "\n+ %i events found ...\n EOF" % (len(coincidences_dict)) + bcolors.ENDC)




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


def runHypo71(coincidences_dict, maxgap=360):
  #pathHypo71 = "%s/hypo71" % (os.environ['PYGEMADIR']) # editar solo en caso de tener claro lo que haces!
  pathHypo71 = "hypo71" 
  phfile = "picks.txt"
  evfile = "eq.pha"

  if os.path.isfile(pathHypo71+"/"+phfile):
    os.remove(pathHypo71+"/"+phfile)

  if os.path.isfile(pathHypo71+"/"+evfile):
    os.remove(pathHypo71+"/"+evfile)

  # generate pickfile from dictionary in hypo71 format
  export_picksfile(coincidences_dict, pickfile=pathHypo71+"/"+phfile)

  # move to hypo71 work folder
  currdir = os.getcwd()
  os.chdir(pathHypo71)





  ######## bassicly is the do_job.py script from Sippl Summer School 2020, University of Concepcion ########

  # read picks file (from "auto-detector")
  picks = hypo71.prepare_picks(phfile)  

  # loop over each event in picks file
  events_dict = {}
  evnum = 1
  for npick in picks:
    try:
      # generate input file "hypo71_input"
      hypo71.generate_input(".", picks[npick]["hypolines"], "info_file", "velmod.hdr") 

      # run hypo71
      hypo71.call_Hypo71(".")

      # read hypo71 results
      out = hypo71.read_output(".")

      # write out in HYPODD format if condition is satisficed
      if float(out["gap"])<=maxgap:
        hypo71.write_pha(out, evfile, evnum)

        # print hypo71 results
        pattern = "[%i] %s  %.4f  %.4f %.1f km    nphases = %i   gap = %.2f  rms = %.4f  0.0 " % (evnum, out["origin_time"].strftime("%Y-%m-%d %H:%M:%S"), float(out["ev_lat"]), float(out["ev_lon"]), float(out["ev_depth"]), int(out["nphases"]), float(out["gap"]), float(out["RMS"])) 
        print(bcolors.BOLD + pattern + bcolors.ENDC)
        for j in out['Tobs']:
          pattern_p = j.replace('_',' ')+' '+out['Tobs'][j]+' 1.00 P'
          #print(pattern_p)
          if j in out['TobsS'].keys():
            pattern_s = j.replace('_',' ')+' '+out['TobsS'][j].strip(' ')+' 1.00 S'
            #print(pattern_s)

        # append data to dictionary
        evtime = out["origin_time"]
        evlat = float(out["ev_lat"])
        evlon = float(out["ev_lon"])
        evdep = float(out["ev_depth"])
        nphases = int(out["nphases"])
        gap = float(out["gap"])
        rms = float(out["RMS"])

        evname = "Event_%03i" % (evnum)
        events_dict.update({evname: {'origin_time'      : evtime, 
                                     'evlat'      : evlat, 
                                     'evlon'      : evlon, 
                                     'evdep'      : evdep, 
                                     'nphases'      : nphases, 
                                     'gap'      : gap, 
                                     'rms'      : rms
                                     } 
                            })

    except:
      continue

    evnum += 1

  # return to your folder
  os.chdir(currdir)

  return events_dict





