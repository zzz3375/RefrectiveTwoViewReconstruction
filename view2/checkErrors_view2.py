#%%
import numpy as np
from DataView2 import point2d as point2ds
from DataView2 import point3d as point3ds
from DataView2 import D0
import matplotlib.pyplot as plt
n = point2ds.shape[0]

#%%
from cameraCalibration import K, R, T
R_T = np.column_stack([R, T])
point_2ds_eular = np.column_stack( [point2ds, np.ones(n)] )
point_3ds_eular = np.column_stack( [point3ds, np.ones(n)] )

#%%
errs = []
# err_rels = []
for i in range(n):
    point2d_eular = point_2ds_eular[i] 
    point3d_eular = point_3ds_eular[i]
    point2d_computed = K @ R_T @ point3d_eular
    point2d_computed = point2d_computed / point2d_computed[-1]
    err_vec = point2d_eular - point2d_computed
    # err_norm = np.linalg.norm(err_vec, 2)
    errs.append( err_vec )
    # point2d_eular_norm = np.linalg.norm(point2d_eular, 2) 
    # err_rels.append( err_norm / point2d_eular_norm )
errs = np.array( errs )
print(errs)
# the error is amazingly under control
# %%
