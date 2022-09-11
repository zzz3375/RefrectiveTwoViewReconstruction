#%%
from cameraCalibration import R, K, T
import numpy as np
from scipy import linalg
# %%
p_2del = np.array([1841, 880, 1])
L1 = R.T @ linalg.inv(K) @ p_2del
L1 = L1 / linalg.norm(L1)
#%%
from cameraCalibration import camera_coordinates as CC
from DataView2 import D0
#%%
Pz = ( 0 -CC[0] )/ L1[0]
test_coordinate = CC + L1 * Pz
actual_coordinate = np.array([0, 8, 6]) * D0
err = test_coordinate - actual_coordinate
print("None refractive err", err )
print("None refractive errNorm", linalg.norm( err ) )
# %%
