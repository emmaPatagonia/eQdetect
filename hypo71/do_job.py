import hypo71
import numpy as np 
import matplotlib.pyplot as plt
import os, glob
from obspy.geodetics.base import degrees2kilometers, calc_vincenty_inverse

# read picks file (from "auto-detector")
pickfile = "picks.txt"
maxgap = 360

picks = hypo71.prepare_picks(pickfile)

# define output file
event_file = "eq.pha"
if os.path.isfile(event_file):
  os.remove(event_file)

# loop over each event in picks file
count = 1
for npick in picks:
  # generate input file "hypo71_input"
  hypo71.generate_input(".", picks[npick]["hypolines"], "info_file", "velmod.hdr") 
  # run hypo71
  hypo71.call_Hypo71(".")
  # read hypo71 results
  out = hypo71.read_output(".")


  # write out in HYPODD format if condition is satisficed
  if float(out["gap"])<=maxgap:
    hypo71.write_pha(out, event_file, count)

    # print hypo71 results
    pattern = "%s  %.4f  %.4f %.1f    %i   %.2f  %.4f  0.0 %i" % (out["origin_time"].strftime("%Y-%m-%d %H:%M:%S"), float(out["ev_lat"]), float(out["ev_lon"]), float(out["ev_depth"]), int(out["nphases"]), float(out["gap"]), float(out["RMS"]), count) 
    print(pattern)
    for j in out['Tobs']:
      print(j.replace('_',' ')+' '+out['Tobs'][j]+' 1.00 P')
      if j in out['TobsS'].keys():
        print(j.replace('_',' ')+' '+out['TobsS'][j].strip(' ')+' 1.00 S')

    count += 1





"""
    # BOOTSTRAP: statistical analysis of perturbation of travel times
    do_bootstrap = False
    if do_bootstrap:
      max_perturbation = 0.2 # in seconds
      num_iters = 1000

      origin_time = out["origin_time"]
      evlat = float(out["ev_lat"])
      evlon = float(out["ev_lon"])
      evdep = float(out["ev_depth"])

      # init loop iterations
      diffs = []
      for i in range(num_iters):
        print("iteration #%i (bootstrap)..."%(i+1))
        # compute new arrival times pertubed
        outfile = open("tmp.txt", "w")
        outfile.write("#"*40)
        outfile.write("\n")
        for stname in out["Tobs"]:
          rnd_p = np.random.normal(loc=0.0, scale=max_perturbation)
          print(rnd_p)
          tp_new = out["origin_time"] + float(out["Tobs"][stname]) + rnd_p
          pattern = "%s P %.2f %s" % (stname, tp_new.timestamp, out['static_weight'][stname])
          if stname in out['TobsS'].keys():
            rnd_s = np.random.normal(loc=0.0, scale=max_perturbation)
            ts_new = out["origin_time"] + float(out["TobsS"][stname]) + rnd_s
            pattern += " S %.2f %s" % (ts_new.timestamp, out['Sweight'][stname])
          pattern += "\n"
          outfile.write(pattern)
        outfile.write("\n")
        outfile.close()

        # compute new hypocenter
        new_picks = hypo71.prepare_picks("tmp.txt")
        hypo71.generate_input(".", new_picks[1]["hypolines"], "info_file", "velmod.hdr") 
        hypo71.call_Hypo71(".")
        new_out = hypo71.read_output(".")
        os.remove("tmp.txt")

        new_origin_time = new_out["origin_time"]
        new_evlat = float(new_out["ev_lat"])
        new_evlon = float(new_out["ev_lon"])
        new_evdep = float(new_out["ev_depth"])

        dx = degrees2kilometers( new_evlon - evlon )
        dy = degrees2kilometers( new_evlat - evlat )
        dz = new_evdep - evdep
        dist = calc_vincenty_inverse(evlat, evlon, new_evlat, new_evlon)[0]/1000
        #std = 
        diffs.append([dx,dy,dz,dist])
      diffs = np.array(diffs)

      # show results
      fig = plt.figure()
      ax = fig.add_subplot(1, 1, 1)
      ax.spines['left'].set_position('center')
      ax.spines['bottom'].set_position('center')
      ax.spines['right'].set_color('none')
      ax.spines['top'].set_color('none')
      ax.xaxis.set_ticks_position('bottom')
      ax.yaxis.set_ticks_position('left')
      ax.minorticks_on()
      ax.set_xlim(-50,50)
      ax.set_ylim(-50,50)

      ax.scatter(diffs.T[0], diffs.T[1], s=25, marker="o", linewidths=0.8, alpha=0.8, c=diffs.T[3])

      plt.show()
"""
