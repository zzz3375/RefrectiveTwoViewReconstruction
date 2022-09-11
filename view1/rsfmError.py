#%%
import numpy as np
from scipy import linalg
from cameraCalibration import T, R, K
from cameraCalibration import camera_coordinates as CC
from DataView1 import D0, Q0x

 
#%%

p_2del = np.array([1507, 830, 1])
L1 = R.T @ linalg.inv( K ) @ p_2del
L1 = L1 / linalg.norm(L1)
t = ( Q0x - CC[0] ) / L1[0]
Q0 = CC + t * L1
#%%

lx, ly, lz = L1
lamta = 1/1.33
lrx = np.sqrt( 1 - lamta**2 + lamta**2 * lx**2) if lx > 0 else -np.sqrt( 1 - lamta**2 + lamta**2 * lx**2)
lry = lamta * ly
lrz = lamta * lz
L1r = np.array([lrx, lry, lrz])
#   refractive line L1r = Q + L1r * t ( t is a variable )


t0 = (0-Q0[0]) / L1r[0]
test_coordinate = Q0 + L1r * t0
actual_coordinate = np.array([0, 2, 6]) * D0 
err = test_coordinate - actual_coordinate
print("one water layer refractive err", err ) 
print("one water layer refractive erNorm", linalg.norm( err ) )
# the error is amazingly under control, but not perfect


 








# %%
