import numpy as np
import cv2

img=np.zeros([2000,2000])
for i in [0,2,4,6,8]:
    for j in [0,2,4,6,8]:
        img[i*200:(i+1)*200,j*200:(j+1)*200]=255

for i in [1,3,5,7,9]:
    for j in [1,3,5,7,9]:
        img[i*200:(i+1)*200,j*200:(j+1)*200]=255

cv2.imshow("win1",img)
cv2.waitKey()
cv2.destroyAllWindows()
cv2.imwrite("chessboard.png",img)
