import sys, re
import h5py
import numpy as np
from scipy import interpolate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

MAXERROR = 0.05

matplotlib.rcParams['font.size'] = 18
matplotlib.rcParams['font.family'] = 'serif'
#matplotlib.rc('text', usetex=True)
matplotlib.rcParams['xtick.major.size'] = 7
matplotlib.rcParams['xtick.major.width'] = 2
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['xtick.minor.size'] = 4
matplotlib.rcParams['xtick.minor.width'] = 2
matplotlib.rcParams['ytick.major.size'] = 7
matplotlib.rcParams['ytick.major.width'] = 2
matplotlib.rcParams['ytick.minor.size'] = 4
matplotlib.rcParams['ytick.minor.width'] = 2
matplotlib.rcParams['axes.linewidth'] = 2


filebasename = sys.argv[1]
nstart   = int(sys.argv[2])
nstop    = int(sys.argv[3])

filenamelist = [filebasename + "%04d"% x for x in range(nstart,nstop+1)]

filenamelist = filenamelist
print(filenamelist)
colorlist = ["red","orange","gold","green","blue","purple"]
labellist = ["Step 5","Step 10","Step 15","Step 20","Step 25","Step 30"]

### Constants
Rsun = 6.96e10
pc   = 3.08e18
clight = 3.0e10

### Data field name
particle_field = 'tracer particles'
real_scalar_field = 'real scalars'
real_runtime_field = 'real runtime parameters'
real_runtime_parameters_field = 'real runtime parameters'
bounding_box_field = 'bounding box'
block_size = 'block size'
node_type_field = 'node type'

### Parameters
ener_thres = 2.0e-15
ener_id = 1 #2
nump_id = 8 #7 #9
posx_id = 9 #8 #10
posy_id = 10 #9 #11
posz_id = 11 #10 #12
velx_id = 15 #14 #16
vely_id = 16 #15 #17
velz_id = 17 #16 #18
nbin = 64
ymax = 3.08e18 #1.62545782181e+13 
xmax = 3.08e18 #1.62545782181e+12 
Rmax = ymax / 2.0
rho_0 = 3.246753247e-17 #1.0e-12
kappa_0 = 1.0 #40.0
flux_0 = 30859982970.0
erad_init = 1.0e10
Dim = 3.0
D_diff = clight / (Dim * kappa_0 * rho_0)

### Defining opacity table 
### Column IDs
logrho_id   = 0
logT_id     = 1
logPgas_id  = 2
logkappa_id = 3

### Reading opacity tables
#infile = np.loadtxt("/home1/02276/bthtsang/py_scripts/kappa_grid.dat", comments="#")

### extracting plot arrays
#k_logrho   = infile[:,logrho_id]
#k_logT     = infile[:,logT_id]
#k_logPgas  = infile[:,logPgas_id]
#k_logkappa = infile[:,logkappa_id]

#k_logrho_grid = np.unique(k_logrho)
#k_logT_grid   = np.unique(k_logT)
#print "kappa table size", k_logrho_grid.size, k_logT_grid.size

#kappa_2d = k_logkappa.reshape((k_logrho_grid.size, k_logT_grid.size), order='C')

### Generating interpolation function
#kappa = interpolate.RectBivariateSpline(k_logrho_grid, k_logT_grid, kappa_2d)
### 
### End opacity table

### creating urad(z) plot
fig = plt.figure()
ax = fig.add_subplot(111)
for i in range(len(filenamelist)):
  filename = filenamelist[i]
  print(filename)
  inputfilenum = re.findall('\d\d\d\d', filename)[0]
  f = h5py.File(filename, 'r')

  reals = f[real_scalar_field]
  real_scalars_dict = {real[0].decode('ascii').strip():real[1] for real in reals}

  time = real_scalars_dict['time']

  real_runs = f[real_runtime_field]
  real_runs_dict = {real[0].decode('ascii').strip():real[1] for real in real_runs}

  xmin = real_runs_dict['xmin']
  xmax = real_runs_dict['xmax']
  ymin = real_runs_dict['ymin']
  ymax = real_runs_dict['ymax']
  zmin = real_runs_dict['zmin']
  zmax = real_runs_dict['zmax']

  domain_center = np.array([0.5*(xmax - xmin), 0.5*(ymax - ymin), 0.5*(zmax - zmin)]) 

  ### Computing cell coordinates
  

  ### Only one type of radiation particles
  mcps = f[particle_field]
  mcps_array = np.array(mcps).T

  mcps_nump = mcps_array[nump_id]
  mcps_ener = mcps_array[ener_id]
  mcps_erad = mcps_nump * mcps_ener

  mcps_posx = mcps_array[posx_id]
  mcps_posy = mcps_array[posy_id]
  mcps_posz = mcps_array[posz_id]

  mcps_velx = mcps_array[velx_id]
  mcps_vely = mcps_array[vely_id]
  mcps_velz = mcps_array[velz_id]

  ### Create plot arrays
  #real_runtime_parameters = f[real_runtime_parameters_field]
  #print real_runtime_parameters[1]
  dR = Rmax / float(nbin)  # sph symmetric, R = 0.5 pc
  R_edges = np.arange(0.0,Rmax+0.1, dR)
  R_in    = R_edges[:-1]
  R_out   = R_edges[1:]
  R_centers = R_edges[:-1] + 0.5 * dR
  R_centers_plot = R_centers / pc
  R_centers_norm = R_centers / Rmax

  volume  = (4./3.)*np.pi*(R_out**3 - R_in**3)

  ### List for urad
  erad_list = np.zeros(R_centers.size)

  ### Lists for mu
  mu_p_list = np.zeros(R_centers.size)
  num_p_list = np.zeros(R_centers.size)

  mu_m_list = np.zeros(R_centers.size)
  num_m_list = np.zeros(R_centers.size)

  ### Lists for mean positions
  posx_list = np.zeros(mcps_erad.size)
  posy_list = np.zeros(mcps_erad.size)
  posz_list = np.zeros(mcps_erad.size)

  ### Collecting MCPs
  for ii in range(mcps_erad.size):
    erad = mcps_erad[ii]

    posx = mcps_posx[ii] - domain_center[0]
    posy = mcps_posy[ii] - domain_center[1]
    posz = mcps_posz[ii] - domain_center[2]
    radius = np.sqrt(posx*posx + posy*posy + posz*posz)
    r_id = int(np.floor(radius / dR))

    if (r_id < nbin-1):
      posx_list[ii] = posx   
      posy_list[ii] = posy
      posz_list[ii] = posz 
 
      erad_list[r_id] += erad

      vely = mcps_vely[ii]
      mu   = vely / clight

      if (mu >= 0.0):
        mu_p_list[r_id] += mu
        num_p_list[r_id] += 1.0
      else:
        mu_m_list[r_id] += mu
        num_m_list[r_id] += 1.0

  ### Report average value
  posx_mean = np.mean(posx_list)
  posy_mean = np.mean(posy_list)
  posz_mean = np.mean(posz_list)
  print ("Mean position", posx_mean, posy_mean, posz_mean)

  ### Lists for urad
  #urad_list = erad_list / (dz * xmax)
  urad_list = erad_list / volume
  urad_list_norm = urad_list * Rmax**3 / erad_init
  urad_ana  = np.zeros(R_centers.size)

  urad_0 = 3.0*kappa_0*rho_0*flux_0*ymax/clight
  #urad_ana = urad_0 - 3.0*kappa_0*rho_0*flux_0*R_centers/clight
  urad_ana = erad_init / (4.0*np.pi*D_diff*time)**(Dim/2.0)
  urad_ana *= np.exp(-R_centers**2/(4.0*D_diff*time))
  urad_ana_norm = urad_ana * Rmax**3 / erad_init

  #avg_mu_p = mu_p_list / num_p_list
  #avg_mu_m = mu_m_list / num_m_list


  #ax.plot(R_centers_plot, urad_list, 'k.')
  #ax.plot(R_centers_plot, urad_ana,  'k-')
  ax.plot(R_centers_norm, urad_ana_norm, color=colorlist[i],linestyle='-',label=labellist[i])
  ax.scatter(R_centers_norm, urad_list_norm, color=colorlist[i])

  # calculate errors
  error = np.sum(np.abs(urad_ana_norm - urad_list_norm)) / np.sum(urad_ana_norm)
  
ax.set_ylim(0,200)
ax.set_xlim(0,0.5)

trash = np.array([-1,-2])
plt.plot(trash,trash,'k-',label="Analytic")
plt.plot(trash,trash,'k.',label="Monte Carlo")
plt.ylabel(r'Radiation Energy Density (erg cm$^{-3}$)')
plt.xlabel(r'$r/R$')
plt.minorticks_on()
plt.legend(frameon=False)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(which='both',direction='in')
            
outfilename = "mcp_urad.pdf"
plt.savefig(outfilename,dpi=300,bbox_inches='tight')

print("error = ",error)
if(error > MAXERROR):
  sys.exit(1)
### Creating average mu plot
#  fig = plt.figure()
#  ax = fig.add_subplot(111)

#  ax.plot(R_centers_plot, avg_mu_p, 'k.')
#  ax.plot(R_centers_plot, avg_mu_m, 'b.')
#  ax.axhline(0.0, linestyle='--')
#  ax.set_ylim(-1.0, 1.0)

#  plt.ylabel(r'Average Directional Cosine $\mu$')
#  plt.xlabel(r'Distance (cm)')

#  outfilename = "mcp_mu_" + inputfilenum + ".pdf"
#  plt.savefig(outfilename,dpi=300,bbox_inches='tight')

