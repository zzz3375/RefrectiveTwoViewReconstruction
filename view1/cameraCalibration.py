#                                   reinventing the wheel
#%%
import numpy as np
from DataView1 import point2d, point3d
from scipy import linalg
#%%                                 transfer to eular coordinates
n = point2d.shape[0]
point3d_eular = np.ones([n,4])
point3d_eular[:,:3] = point3d


#%%                                 create P matrix
P = np.zeros([2*n,12])

for i in range(n):
    ui, vi = point2d[i]
    Pi = point3d_eular[i].reshape([4,1])
    P[2*i, :4 ] = Pi.T.ravel()
    P[2*i, -4:] = -ui*Pi.T.ravel()
    P[2*i+1, 4:8] = Pi.T.ravel()
    P[2*i+1, -4:] = -vi*Pi.T.ravel()



#%%                                   single value decomposition
U, s, Vh = linalg.svd(P)

#%%                                      choose Vh[-1] as m
m = Vh[-1].reshape((3,4))

a1,a2,a3 = m[:,:3]
ro2 = 1/(a3 @ a3)
ro = np.sqrt(ro2)

u0 = (a1 @ a3 )/(a3 @ a3) 
v0 = (a2 @ a3 )/(a3 @ a3)

cos_theta = (( np.cross(a1,a3) @ np.cross(a2,a3) ) / np.linalg.norm(np.cross(a1,a3), 2)/np.linalg.norm(np.cross(a2,a3), 2))
theta = np.arccos( cos_theta )

alpha = ro2*np.linalg.norm(np.cross(a1,a3), 2)*np.sin(theta)
beta = ro2*np.linalg.norm(np.cross(a2,a3), 2)*np.sin(theta)

r1 =  -np.cross(a2, a3) / np.linalg.norm(np.cross(a2, a3))
r3 =  ro * a3
r2 =  -np.cross(r3,r1)
#%%
K = np.array([[alpha, -alpha/np.tan(theta), u0],
              [0, beta/np.sin(theta), v0],
              [0,             0,        1]])
T = ro * np.linalg.inv(K) @ m[:,-1]
R = np.array([r1,r2,r3])
camera_coordinates = linalg.inv(R) @ (-T)
# %%
np.save("K1.npy", K)
np.save("T1.npy", T)
np.save("R1.npy", R)