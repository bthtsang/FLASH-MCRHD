import numpy as np
import h5py
import glob
import matplotlib.pyplot as plt

stem = "cartimc_hdf5_chk_"
filenames = sorted(glob.glob(stem+"*"))
a = 7.5657e-15 # erg cm^-3 K^-1

def plot_minmax(name):
    minval = []
    maxval = []
    time = []
    for filename in filenames:
        f = h5py.File(filename,"r")

        time.append( f["real scalars"][0][1] )
    
        val = np.array(f[name])
        if name=="urad":
            val = (val/a)**.25 # radiation temperature
        minval.append(val.min())
        maxval.append(val.max())

        f.close()

    iter = np.arange(len(time))
    plt.clf()
    plt.plot(iter,minval)
    plt.plot(iter,maxval)
    plt.savefig("cartradeqm_"+name+".pdf")

def plot_bothtemperatures():
    minval_fluid = []
    maxval_fluid = []
    minval_rad = []
    maxval_rad = []
    time = []
    for filename in filenames:
        f = h5py.File(filename,"r")

        time.append( f["real scalars"][0][1] )
    
        Tfluid = np.array(f["temp"])
        Trad = (np.array(f["urad"])/a)**.25
        minval_fluid.append(Tfluid.min())
        maxval_fluid.append(Tfluid.max())
        minval_rad.append(Trad.min())
        maxval_rad.append(Trad.max())

        f.close()

    iter = np.arange(len(time))
    plt.clf()
    plt.xlabel("t (s)")
    plt.ylabel("T (K)")
    plt.plot(iter,minval_fluid)
    plt.plot(iter,maxval_fluid)
    plt.plot(iter,minval_rad)
    plt.plot(iter,maxval_rad)
    plt.savefig("cartradeqm_TfluidTrad.pdf", bbox_inches="tight")

    # test final temperature
    radSum = np.sum(Trad)
    fluidSum = np.sum(Tfluid)
    error = np.abs((radSum - fluidSum) / (radSum+fluidSum))

    
plot_bothtemperatures()
#plot_minmax("temp")
#plot_minmax("urad")
#plot_minmax("ye  ")
#plot_minmax("flec")
#plot_minmax("eint")

#nparticles = []
#time = []
#for filename in filenames:
#    f = h5py.File(filename,"r")
#
#    time.append( f["real scalars"][0][1] )
#
#    if "tracer particles" in f:
#        nparticles.append( np.shape(f["tracer particles"])[0] )
#    else:
#        nparticles.append(0)
#    f.close()
#
#iter = np.arange(len(time))
#
#plt.clf()
#plt.plot(iter, nparticles)
#plt.savefig("cartradeqm_nparticles.pdf")
