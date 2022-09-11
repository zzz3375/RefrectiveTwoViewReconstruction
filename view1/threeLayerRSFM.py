import numpy as np
from scipy import linalg
from cameraCalibration import T, R, K
from cameraCalibration import camera_coordinates as CC
from DataView1 import D0, Q0x, Q1x, Q2x, Q3x
from DataView1 import n_air, n_glass, n_water

p_2del = np.array([1507, 830, 1])
L0 = R.T @ linalg.inv( K ) @ p_2del
L0 = L0 / linalg.norm(L0)
Pz = ( Q0x - CC[0] ) / L0[0]
Q0 = CC + Pz * L0



lx, ly, lz = L0
lamta = n_air / n_glass
lrx = -np.sqrt( 1 - lamta**2 + lamta**2 * lx**2)
lry = lamta * ly
lrz = lamta * lz
L1 = np.array([lrx, lry, lrz])
# L1 ~ Q0 + L1*t   |L1| = 1



t = (Q1x-Q0x) / L1[0]
Q1 = L1 * t + Q0
lamta = n_glass / n_water
lx, ly, lz = L1
lrx = -np.sqrt( 1 - lamta**2 + lamta**2 * lx**2)
lry = lamta * ly
lrz = lamta * lz
L2 = np.array([lrx, lry, lrz])


t = ( Q2x - Q1x ) / L2[0]
Q2 = Q1 + L2 *t
lamta = n_water / n_glass
lx, ly, lz = L2
lrx = -np.sqrt( 1 - lamta**2 + lamta**2 * lx**2)
lry = lamta * ly
lrz = lamta * lz
L3 = L2 = np.array([lrx, lry, lrz])


t = ( Q3x - Q2x ) / L3[0]
Q3 = Q2 + t * L3
test_coordinate = Q3
actual_coordinate = np.array([0, 2, 6]) * D0 
err = test_coordinate - actual_coordinate
print("multilayer refractive err", err )  


## multilayer can not see an obivious improve