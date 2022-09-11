#%%
import numpy as np
from scipy import linalg
from PointsForReconstruction import point3ds, point2ds_view1, point2ds_view2
from view1.DataView1 import Q0x
import matplotlib.pyplot as plt

K1 = np.load("view1\K1.npy")
R1 = np.load("view1/R1.npy")
T1 = np.load("view1/T1.npy")
K2 = np.load("view2/K2.npy")
R2 = np.load("view2/R2.npy")
T2 = np.load("view2/T2.npy")
#%%
if point2ds_view1.shape[0] != point2ds_view2.shape[0]:
    print("image cannot fully match ! EXITING NOW")
    exit()
else: print("image numbers match well")

n = point2ds_view1.shape[0]
point2ds_el_view1 = np.column_stack([ point2ds_view1, np.ones(n)])    
point2ds_el_view2 = np.column_stack([ point2ds_view2, np.ones(n)])

CC1 = R1.T @ (-T1)
CC2 = R2.T @ (-T2)
#%%
errs_r = np.zeros([n, 3])
errs = np.zeros([n, 3])
alpha1s = np.zeros(n)

for i in range(n):
    # view1 ray spread in air
    p1_2del = point2ds_el_view1[i]
    L1 = R1.T @ linalg.inv(K1) @ p1_2del
    L1 = L1 / linalg.norm(L1)
    
    t1 = ( Q0x - CC1[0] ) / L1[0]
    Q01 = L1 * t1 +CC1
    # view1 refraction
    lamta = 1 / 1.33
    lx1, ly1, lz1= L1
    lrx1_norm = np.sqrt( 1 - lamta**2 + lamta**2 * lx1**2) 
    lrx1 = lrx1_norm if lx1 > 0 else -lrx1_norm
    lry1 = lamta * ly1
    lrz1 = lamta * lz1
    L1r = np.array([lrx1, lry1, lrz1])
    
    # view2 ray spread in air
    p2_2del = point2ds_el_view2[i]
    L2 = R2.T @ linalg.inv(K2) @ p2_2del
    L2 = L2 / linalg.norm(L2)
    
    t2 = ( Q0x - CC2[0] ) / L2[0]
    Q02 = L2 * t2 +CC2
    # view2 refraction
    lamta = 1/1.33
    lx2, ly2, lz2 = L2
    lrx2_norm = np.sqrt( 1 - lamta**2 + lamta**2 * lx2**2)
    lrx2 = lrx2_norm if lx2 > 0 else - lrx2_norm
    lry2 = lamta * ly2
    lrz2 = lamta * lz2
    L2r = np.array([lrx2, lry2, lrz2])
    
    Q12 = Q01 - Q02
    E1 = np.cross(Q12, L2r)
    F1 = np.cross(L1r, L2r)
    t1 = - ( F1 @ E1 ) / ( F1 @ F1 )
    d1 = linalg.norm(E1 + F1 * t1) # min distance of 2 rays
    J1 = Q01 + L1r * t1 # min distance point in ray 1
    # get correct version of L1r and L2r

    Q21 = Q02 - Q01
    E2 = np.cross(Q21, L1r)
    F2 = np.cross(L2r, L1r)
    t2 = - ( F2 @ E2 ) / ( F2 @ F2 )
    d2 = linalg.norm(E2 + F2 * t2)
    J2 = Q02 + L2r * t2

    x = (J1 + J2) / 2
    actual_coordinate = point3ds[i]
    err_r = x - actual_coordinate
    errs_r [i] = err_r
    
    # direct sfm method
    # view1
    CC12 = CC1 - CC2
    E1_d = np.cross(CC12, L2)
    F1_d = np.cross(L1, L2)
    t1_d = - ( F1_d @ E1_d ) / ( F1_d @ F1_d )
    d1_d = linalg.norm(E1_d + F1_d * t1_d)
    J1_d = CC1 + L1 * t1_d
    #view2
    CC21 = CC2 - CC1
    E2_d = np.cross(CC21, L1)
    F2_d = np.cross(L2,L1)
    t2_d = - ( F2_d @ E2_d ) / ( F2_d @ F2_d )
    d2_d = linalg.norm(E2_d + F2_d * t2_d)
    J2_d = CC2 + L2 * t2_d
    x_d = (J1_d + J2_d)/2
    err = x_d - actual_coordinate
    errs[i] = err
    O1 = CC1 - point3ds[i]
    O2 = CC2 - point3ds[i]
    cos_alpha1 = (O1 @ O2) / linalg.norm(O1) / linalg.norm(O2)
    alpha1 = np.arccos(cos_alpha1)
    alpha1 = alpha1/np.pi *180
    alpha1s[i] =alpha1


#%%

plt.scatter(alpha1s, linalg.norm(errs, axis = 1),label = "sfm")
plt.scatter(alpha1s, linalg.norm(errs_r, axis = 1),label = "rsfm")
plt.legend()
plt.xlabel("parallax angle (degree)")
plt.ylabel("translation err (mm)")
plt.savefig("TranslationError.svg",format='svg')


# %%
