
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.size'] = 15


# In[2]:


from qprop.core import Qprop20

q = Qprop20('.')

N_p = q.propagateParam['num-of-particle']
N_r_dim = 3

import numpy as np

traj_r_t_arr = np.fromfile('traj-r-t.bin')
traj_v_t_arr = np.fromfile('traj-v-t.bin')

N_t, _res = divmod(traj_r_t_arr.size, N_p * N_r_dim)
if (_res != 0):
    print("num_of_entry: {} / N_t: {} / _res: {}".format(traj_r_t_arr.size, N_t, _res))


from tdse.coordinate import spherical_to_cartesian
traj_r_t_3d_arr = traj_r_t_arr.reshape(N_t,N_r_dim,N_p)
traj_v_t_3d_arr = traj_v_t_arr.reshape(N_t,N_r_dim,N_p)
x_r_t_3d_arr, y_r_t_3d_arr, z_r_t_3d_arr = spherical_to_cartesian(traj_r_t_3d_arr[:,0,:], traj_r_t_3d_arr[:,1,:], traj_r_t_3d_arr[:,2,:])
x_v_t_3d_arr, y_v_t_3d_arr, z_v_t_3d_arr = spherical_to_cartesian(traj_v_t_3d_arr[:,0,:], traj_v_t_3d_arr[:,1,:], traj_v_t_3d_arr[:,2,:])


# In[3]:


r_max = max([np.abs(x_r_t_3d_arr).max(), np.abs(y_r_t_3d_arr).max(), np.abs(z_r_t_3d_arr).max()]) * 1.1
v_max = max([np.abs(x_v_t_3d_arr).max(), np.abs(y_v_t_3d_arr).max(), np.abs(z_v_t_3d_arr).max()]) * 1.1


# In[4]:

ncol = 2
nrow = 2


fig, axes = plt.subplots(nrows=nrow,ncols=ncol,figsize=(6*ncol,7*nrow))


ax = axes[0,0]
for i in range(N_p):
    line, = ax.plot(x_r_t_3d_arr[:,i], z_r_t_3d_arr[:,i])
    ax.plot(x_r_t_3d_arr[0,i], z_r_t_3d_arr[0,i], '.', color=line.get_color())
    ax.plot(x_r_t_3d_arr[-1,i], z_r_t_3d_arr[-1,i], 'x', color=line.get_color())

ax.axis('square')
ax.set_xlim(-r_max,r_max)
ax.set_ylim(-r_max,r_max)
ax.set_xlabel(r"$x$ / a.u.")
ax.set_ylabel(r"$z$ / a.u.")
ax.set_title("trajectory @ xz-plane")


ax = axes[0,1]
for i in range(N_p):
    line, = ax.plot(x_r_t_3d_arr[:,i], y_r_t_3d_arr[:,i])
    ax.plot(x_r_t_3d_arr[0,i], y_r_t_3d_arr[0,i], '.', color=line.get_color())
    ax.plot(x_r_t_3d_arr[-1,i], y_r_t_3d_arr[-1,i], 'x', color=line.get_color())

ax.axis('square')
ax.set_xlim(-r_max,r_max)
ax.set_ylim(-r_max,r_max)
ax.set_xlabel(r"$x$ / a.u.")
ax.set_ylabel(r"$y$ / a.u.")
ax.set_title("trajectory @ xy-plane")


ax = axes[1,0]
for i in range(N_p):
    line_v, = ax.plot(x_v_t_3d_arr[:,i], z_v_t_3d_arr[:,i])
    ax.plot(x_v_t_3d_arr[0,i], z_v_t_3d_arr[0,i], '.', color=line_v.get_color())
    ax.plot(x_v_t_3d_arr[-1,i], z_v_t_3d_arr[-1,i], 'x', color=line_v.get_color())

ax.axis('square')
ax.set_xlim(-v_max,v_max)
ax.set_ylim(-v_max,v_max)
ax.set_xlabel(r"$v_x$ / a.u.")
ax.set_ylabel(r"$v_z$ / a.u.")
ax.set_title("velocity @ xz-plane")


ax = axes[1,1]
for i in range(N_p):
    line_v, = ax.plot(x_v_t_3d_arr[:,i], y_v_t_3d_arr[:,i])
    ax.plot(x_v_t_3d_arr[0,i], y_v_t_3d_arr[0,i], '.', color=line_v.get_color())
    ax.plot(x_v_t_3d_arr[-1,i], y_v_t_3d_arr[-1,i], 'x', color=line_v.get_color())

ax.axis('square')
ax.set_xlim(-v_max,v_max)
ax.set_ylim(-v_max,v_max)
ax.set_xlabel(r"$v_x$ / a.u.")
ax.set_ylabel(r"$v_y$ / a.u.")
ax.set_title("velocity @ xy-plane")
    
fig.tight_layout()


# In[6]:


fig.savefig("traj-r_t-and-v_t.png")

