#module hypo71
#coding=utf-8

from obspy.core import UTCDateTime
import subprocess, os
from subprocess import Popen
from pylab import *

"""
module for calling HYPO71 from Python and evaluating its output
call function find_location for location of a single event (see calling script hypo71_call.py)
if necessary, iterations over different starting depths and the elimination of the picks with the biggest residuals are performed
"""

def get_stat_values(info_file):
  """
  function that reads station names and coordinates from info_file and turns them into the format needed for the HYPO71 input file
  info_file: path to the file containing the station information
  returns a list of strings readymade for the HYPO71.inp file
  """
  cmd = 'awk \'/^/ {print $1,$8,$9,$10}\' '+info_file
  process = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
  output = process.communicate()[0]
  outlist = output.split(b'\n')
  strlist = []
  for i in range(len(outlist)-1):
    if outlist[i][0] != '#':
      stat_name,stat_lat,stat_lon,stat_ev = outlist[i].split(b' ')
      lat_deg,lat_dec = stat_lat.split(b'.')
      lon_deg,lon_dec = stat_lon.split(b'.')
      #convert decimal degrees into minutes
      minutes_lat = '%05.2f' % (float(b'0.'+lat_dec)*60.)
      minutes_lon = '%05.2f' % (float(b'0.'+lon_dec)*60.)
      #expand station names shorter than 4 characters with _ or __
      stat_name = str(stat_name, 'utf-8')
      if len(stat_name) == 2:
        stat_name = stat_name+'__'
      elif len(stat_name) == 3:
        stat_name = stat_name+'_'
      #determine hemisphere(s)
      if float(lat_deg) >= 0.:
        lat_mark = 'N'
      else:
        lat_mark = 'S'
      if float(lon_deg) >= 0.:
        lon_mark = 'E'
      else:
        lon_mark = 'W'
      #format decimal degrees
      lat_d = '%02i' % (abs(int(lat_deg)))
      lon_d = '%03i' % (abs(int(lon_deg)))

      #format station elevations
      
      stat_e = '%4s' % (str(int(float(stat_ev))))
      if stat_ev == 'XXX':
        stat_e = '   0'
      
      #all into a string for this one station
      stg = '  '+stat_name+lat_d+minutes_lat+lat_mark+lon_d+minutes_lon+lon_mark+stat_e+' 00.00\n'
      #append that string to the list
      strlist.append(stg)
  return strlist

def prepare_picks(phase_file): #indict: Anpassung gemacht!!
  """
  convert phase arrivals into the format read by Hypo71
  phase_file: Path to the file containing the phase picks
  returns dictionaries hypo_dict and nr_dict
  hypo_dict: for each event number, a list of the input strings needed for hypo71.inp (in alphabetical order)
  orig_dict: event header for each event number (for later evaluation)
  """
  phases = open(phase_file,'r')
  ph_data = phases.readlines()
  phases.close()

  indict = {}
  num = 0
  for a in range(len(ph_data)):
    if ph_data[a][0] == '#':
      str_list = []
      indict[num+1] = {}
      indict[num+1]['hypolines'] = []
      for b in range((a+1),(a+2000)):
        #if b < len(ph_data):
          if ph_data[b] == '\n':
            num += 1
            break
          if ph_data[b][0] == '#':
            num += 1
            break
          sflag = False
          try: 
            stat,typP,ppick,weight = ph_data[b].split(None)
          except:
            stat,typP,ppick,weight,typS,spick,weightS = ph_data[b].split(None)
            sflag = True
          picktime = UTCDateTime(float(ppick))
          if sflag:
            picktimeS = UTCDateTime(float(spick))
          if len(stat) == 2:
            stat = stat+'__'
          elif len(stat) == 3:
            stat = stat+'_'
          hyp_year = str(picktime.year)[2:]
          hyp_month = '%02d' % (picktime.month)
          hyp_day = '%02d' % (picktime.day)
          hyp_hour = '%02d' % (picktime.hour)
          hyp_min = '%02d' % (picktime.minute)
          sc = picktime.second
          msec = picktime.microsecond
          t_sec = float(sc + msec/1e6)
          hyp_sec = '%5.2f' % (t_sec)
          if sflag:
            sptime = '%6.2f' % ( (picktimeS - picktime) + picktime.second )
            strg = stat+'IP_'+weight+' '+hyp_year+hyp_month+hyp_day+hyp_hour+hyp_min+hyp_sec+'      '+sptime+'IS_'+weightS+'\n'
          else:
            strg = stat+'IP_'+weight+' '+hyp_year+hyp_month+hyp_day+hyp_hour+hyp_min+hyp_sec+'\n' #weight 0 for all picks
          indict[num+1]['hypolines'].append(strg)
      indict[num]['hypolines'].sort()
  return indict


def prepare_vmod(velmod):
  """
  function to read in and format velmod.hdr for use in hypo71.inp (no big preparation necessary, just reading it in and handing it back...)
  velmod: path to velocity model file (only P velocities, S is calculated via vp/vs-value given in the control cards)
  """
  infile = open(velmod,'r')
  data = infile.readlines()
  infile.close()
  return data

def generate_input(path,phase_dict,info_file,velmod,fdep=35,w1=500,w2=999,vpvs=1.73):
  """
  generates hypo71 input file hypo71.inp
  Control card values to be given as integers (exception: vp/vs-ratio)
  writes file hypo71.inp to given path
  out_path: specifies the path to where the output file hypo71.inp is written
  phase_dict: one entry of hypo_dict (as returned by prepare_picks), i.e. picks for one single event in specific hypo71 input format
  info_file: path to file containing station information
  velmod: location of velmod.hdr file
  fdep: control card entry for starting depth for HYPO71 (default: 35km)
  w1: control card entry for epicentral distance below which stations are fully weighted (default: 200km)
  w2: control card entry for epicentral distance above which stations get 0 distance weight (default: 400km)
  vpvs: vp/vs ratio for calculation of S-velocity model (default: 1.73)
  """
  outfile = open(path+'/hypo71.inp','w')
  # setting the TEST variables
  strg = 'HEAD                     HYPO71PC FOR AUTO GIANT\nRESET TEST(01)=0.100000\nRESET TEST(02)=10.000000\nRESET TEST(03)=2.000000\nRESET TEST(04)=0.050000\nRESET TEST(05)=5.000000\nRESET TEST(06)=4.000000\nRESET TEST(07)=-0.870000\nRESET TEST(08)=2.000000\nRESET TEST(09)=0.003500\nRESET TEST(10)=100.000000\nRESET TEST(11)=20.00000\nRESET TEST(12)=0.500000\nRESET TEST(13)=1.000000\n\n'
  outfile.write(strg)
  #fetching station coordinates and elevations
  statlist = get_stat_values(info_file)
  for g in statlist:
    outfile.write(g)
  #add empty line
  outfile.write('\n')
  #fetching velocity model
  vmod = prepare_vmod(velmod)
  for h in vmod:
    outfile.write(h)
  #another empty line
  outfile.write('\n')
  #setting control cards
  fdep = '%4d' % (fdep)
  vpvs = '%4.2f' %(vpvs)
  controlcards = fdep+'. '+str(w1)+'. '+str(w2)+'. '+str(vpvs)+'    2    1   18         1         1   11                    \n'
  outfile.write(controlcards)
  #fetching phase pick data (get via event number from dictionary)
  for i in phase_dict:
    outfile.write(i)
  #last line
  outfile.write('                 10                                                             \n')
  outfile.close()

def initialize(path,phase_file): #Anpassung indict gemacht !!
  """
  creates piping file and dictionary of phases
  file output: hypo71_input
  hands back dictionaries hypo_dict and orig_dict (definition see function prepare_picks)
  """
  if mode == 'pre':
    indict = prepare_picks(phase_file)
  elif mode =='post' or mode == 'postS':
    indict = read_MPXout(phase_file)
  outfile = open(path+'/hypo71_input','w')
  outfile.write('hypo71.inp\nhypo71.prt\nhypo71.pun\n')
  outfile.close()
  return indict
  

def call_Hypo71(path):
  """
  calls HYPO71 via system call; output written to file hypo71.prt
  """
  inp_file = path+'/hypo71_input'
  cmd = 'hypo71pc < '+inp_file + '> out.log'
  os.system(cmd)

def read_output(path):
  """
  reads hypo71 output from the hypo71.prt file and extracts vital information into a dictionary, which is returned
  prt_file: path to .prt-file to be harvested
  returns loc_dict, a dictionary containing all important information on the obtained location
  """
  prt_file = path+'/hypo71.prt'
  infile = open(prt_file,'r')
  dat = infile.readlines()
  infile.close()
  
  cmd = 'awk \'/^  DATE/ {print NR}\' '+prt_file
  process = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
  output = process.communicate()[0]
  try:
    num = int(output)
    header_line = dat[num-1] #extract hemisphere information
    dum,dum,dum,n,dum,e,dum = header_line.split(None,6)
    sum_line = dat[num]
    #extraction of vital information from summary hypocenter line
    yr = sum_line[1:3]
    if int(yr) < 10:
      year = '200'+yr.strip(' ')
    else:
      year = '20'+yr
    mn = sum_line[3:5]
    day = sum_line[5:7]
    hr = sum_line[8:10]
    min = sum_line[10:12]
    sec = sum_line[13:18]
    origin_time = UTCDateTime(int(year),int(mn),int(day),int(hr),int(min),float(sec))
    # revision needed here -- currently would lead to false results for negative coordinate (i.e. S and W hemispheres)
    lat_deg,lat_min = sum_line[19:27].split('-')
    ev_lat = str(float(lat_deg) + (float(lat_min)/60.))
    lon_deg,lon_min = sum_line[28:37].split('-')
    ev_lon = str(float(lon_deg) + (float(lon_min)/60.))
    ev_depth = sum_line[38:44].strip(' ') #to get rid of leading spaces...
    Nobs = sum_line[52:54].strip(' ')
    gap = sum_line[58:61].strip(' ')
    RMS = sum_line[63:68].strip(' ')
    loc_dict = {}
    if n=='S': ev_lat='-'+ev_lat #adjustment for western and southern hemispheres
    if e=='W': ev_lon='-'+ev_lon
    loc_dict['origin_time'] = origin_time
    loc_dict['ev_lat'] = ev_lat 
    loc_dict['ev_lon'] = ev_lon
    loc_dict['ev_depth'] = ev_depth
    loc_dict['nphases'] = Nobs
    loc_dict['gap'] = gap
  #  loc_dict['RMS'] = RMS  #disabled --> RMS will be determined further down directly from station residuals
    loc_dict['goodness'] = True
    loc_dict['Tobs'] = {}
    loc_dict['TobsS'] = {}
    loc_dict['Sweight'] = {}
    loc_dict['Sres'] = {}
    loc_dict['Sstatw'] = {}
    loc_dict['static_weight'] = {}
    loc_dict['polarity'] = {}  
    loc_dict['P_weight'] = {}
    loc_dict['P_station_res'] = {}
    #getting station residuals from further down in prt file
    res_lst = []
    for i in range((num+3),(num+1000)):
      if dat[i] == '\n' or dat[i][3:7] == 'DATE' or dat[i][0] == '1':
        break
      sta = dat[i][1:5]
      stat_w = dat[i][23:24]
      pol = dat[i][22:23]
      tpobs = dat[i][36:41]
      if tpobs[0] == '0': #should have been three digit...otherwise this field would stay blank...
        tpobs = '1'+tpobs
      pweight = dat[i][61:65]
      res = dat[i][53:59]
      spobs = dat[i][109:115] #become spaces if no S pick exists 
      s_swgt = dat[i][102:103]
      sweight = dat[i][124:127]
      sres = dat[i][115:121]
      loc_dict['P_station_res'][sta] = res
      loc_dict['Tobs'][sta] = tpobs
      loc_dict['static_weight'][sta] = stat_w
      loc_dict['polarity'][sta] = pol
      loc_dict['P_weight'][sta] = pweight
      if spobs.strip(' ') != '': #existence of S phase
        loc_dict['TobsS'][sta] = spobs #relative to origin time
        loc_dict['Sweight'][sta] = sweight
        loc_dict['Sres'][sta] = sres
        loc_dict['Sstatw'][sta] = s_swgt
      try: 
        res_lst.append(float(res))
        res_lst.append(float(sres))
      except:
        pass  
      '''if res == '******':
        del loc_dict['P_weight'][sta]
        del loc_dict['polarity'][sta]
        del loc_dict['static_weight'][sta]
        del loc_dict['Tobs'][sta]
        del loc_dict['P_station_res'][sta]
      '''  #later: should also incorporate S...
    val = 0
    for t in range(len(res_lst)): 
      v = res_lst[t]**2
      val += v
    val /= float(len(res_lst))
    rms = sqrt(val)
    loc_dict['RMS'] = str('%5.2f'%(rms)) #calculation "by hand" because in some cases HYPO71 sets internal weighting for stations with huge residuals to 0 --> those are then not included in determination of HYPO71's RMS
    #print loc_dict
    return loc_dict
  except ValueError: #i.e. no location was obtained, line in hypo71.prt does not exist
    loc_dict = {}
    loc_dict['RMS'] = 99.#RMS dummy
    loc_dict['ev_depth'] = 999. #depth dummy
    loc_dict['goodness'] = False
    return loc_dict

def write_pha(out, event_file, n):
  '''
  write into .pha file for input to HypoDD/ph2dt
  '''
  import os
  if not os.path.isfile(event_file):
    outfile = open(event_file,'w')
  else:
    outfile = open(event_file,'a')
  #n = 1
  #Event line
  outfile.write('# %4i %02i %02i %02i %02i %05.2f %8.4f %8.4f %5.2f  %i %.2f %.4f 0.0' % (out['origin_time'].year,out['origin_time'].month,out['origin_time'].day,out['origin_time'].hour,out['origin_time'].minute,(out['origin_time'].second + (out['origin_time'].microsecond/1e6)),float(out['ev_lat']),float(out['ev_lon']),float(out['ev_depth']), int(out["nphases"]), float(out["gap"]), float(out["RMS"]) ) )
  outfile.write(' '+str(n)+'\n') 
  #Pick lines
  for j in out['Tobs']:
    outfile.write(j.replace('_',' ')+' '+out['Tobs'][j]+' 1.00 P\n')
    if j in out['TobsS'].keys():
      outfile.write(j.replace('_',' ')+' '+out['TobsS'][j].strip(' ')+' 1.00 S\n')

  outfile.close()