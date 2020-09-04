import numpy as np
import matplotlib.pyplot as plt

info = "eq.pha"

out=[]
f = open(info,"r")
for line in f:
  if line[0]=="#":
    evlat = float(line.split(" ")[7])
    evlon = float(line.split(" ")[8])
    evdep = float(line.split(" ")[9])
    out.append([evlat,evlon])
f.close()
out = np.array(out)

stlats, stlons = np.loadtxt("info_file",usecols=(7,8)).T

plt.scatter(out.T[1], out.T[0], marker='o')
plt.scatter(stlons,stlats,marker="v")
plt.show()