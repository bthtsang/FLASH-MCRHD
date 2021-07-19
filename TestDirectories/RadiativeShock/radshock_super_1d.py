#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import h5py
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


# In[2]:


#pwd


# In[6]:


### checkpoint data field names
dens_field = 'dens'
bndb_field = 'bounding box'
pres_field = 'pres'
temp_field = 'temp'
urad_field = 'urad'
abse_field = 'abse'
emie_field = 'emie'
node_field = 'node type'

### Load up last checkpoint file
blk_id  = 0
ener_id = 1
nump_id = 8
posx_id = 9
posy_id = 10
posz_id = 11
velx_id = 15
vely_id = 16
velz_id = 17

Lmin, Lmax = 0.0, 7.0e10
dim = 256
v0 = -20.0e5

dx = (Lmax - Lmin) / dim
dy = (Lmax - Lmin) / dim
dz = (Lmax - Lmin) / dim
nbx, nby, nbz = 8, 8, 8  # FLASH's block size

arad = 7.565723351185839E-015

### Reading in checkpoint file for energy densities

hdf5_list = ["radshock_hdf5_chk_0005", "radshock_hdf5_chk_0010", "radshock_hdf5_chk_0020"]
#hdf5_list = ["radshock_hdf5_chk_0010"]
nfile = len(hdf5_list)

### Initialize empty dens, temp, urad
dens1D = np.zeros((nfile, dim))
temp1D = np.zeros((nfile, dim))
urad1D = np.zeros((nfile, dim))
abse1D = np.zeros((nfile, dim))
emie1D = np.zeros((nfile, dim))

for fid, fname in enumerate(hdf5_list):
    f = h5py.File(fname, 'r')

    ### Current time
    reals = f["real scalars"]
    real_scalars_dict = {real[0].decode('ascii').strip():real[1] for real in reals}
    time = real_scalars_dict['time']
    print ("time", time)

    ### Grid data
    node = np.array(f[node_field])
    dens_data = np.array(f[dens_field])
    temp_data = np.array(f[temp_field])
    urad_data = np.array(f[urad_field])
    abse_data = np.array(f[abse_field])
    emie_data = np.array(f[emie_field])
    bndb_data = np.array(f[bndb_field])
    ### screen only for leaf blocks
    leaf_node = (node == 1)  # 1 means leaf blocks
    num_leaf  = np.count_nonzero(leaf_node)
    print (f"Number of leaf blocks = {num_leaf}")

    dens = dens_data[leaf_node]
    temp = temp_data[leaf_node]
    urad = urad_data[leaf_node]
    abse = abse_data[leaf_node]
    emie = emie_data[leaf_node]
    bndb = bndb_data[leaf_node]

    # loop over blocks to populate the data cubes
    for b in range(num_leaf):
        # inverted xyz indexing because of fortran style HDF5 output
        ix = int((bndb[b][0][0] - Lmin) / dx)  # [0] for lower bound

        # populate one block's data
        dens1D[fid, ix:ix+nbx] = dens[b]
        temp1D[fid, ix:ix+nbx] = temp[b]
        urad1D[fid, ix:ix+nbx] = urad[b]
        abse1D[fid, ix:ix+nbx] = abse[b]
        emie1D[fid, ix:ix+nbx] = emie[b]

# print (temp1D)
    
trad1D = (urad1D/arad)**(0.25)
print ("temp", np.amin(temp1D), np.amax(temp1D))
print ("trad", np.amin(trad1D), np.amax(trad1D))
    
### MCP data
mcps = f['tracer particles']
mcp_ener = mcps[:, ener_id]
mcp_nump = mcps[:, nump_id]
mcp_ener_eV = mcp_ener



### Closing the checkpoint file
f.close()


# In[7]:


### Make the temperature vs x plot

### Preparing the average array, average over y and z directions

cell_count = np.linspace(1,dim,dim)-0.5
posx = cell_count * dx #- v0 * time
posx = posx/1.0e10


fig, ax = plt.subplots(figsize=(8,6))

color_list = ["C0", "C1", "C2"]
### Plot each time step
for fid in range(nfile):
    ax.plot(posx, temp1D[fid,:], f"{color_list[fid]}-", label="Gas")
    ax.plot(posx, trad1D[fid,:], f"{color_list[fid]}--", label="Rad.")

# ax.plot(posx, abse1D, "C0-", label="Abs.")
# ax.plot(posx, np.abs(emie1D), "C1-", label="Emi.")
ax.set_ylabel('Temperature (K)', fontsize=28)
ax.set_xlabel(r'x (10$^{10}$ cm)', fontsize=28)

#ax.set_xlim(2.0, 5.0)
ax.set_xlim(0.0, 7.0)

ax.tick_params(axis='both', direction='in', which='both', labelsize=20, right=True, top=True)

# ax.xaxis.set_major_locator(MultipleLocator(2))
# ax.xaxis.set_minor_locator(MultipleLocator(1))
# ax.yaxis.set_major_locator(MultipleLocator(250))
# ax.yaxis.set_minor_locator(MultipleLocator(50))

plt.tick_params(axis='both', labelsize=24)
plt.legend(fontsize=18, loc=1, frameon=False, ncol=1)

plt.tight_layout()
plt.savefig(f"temp_vs_x.pdf")


