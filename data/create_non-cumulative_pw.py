import numpy as np

data = np.load("DATA-2015-12-31-00.npz")
pw = data["pw"]
right_pw = np.zeros((49,180,360))
for i in range(1,47):
    for j in range(0,179):
        for k in range(0,359):
            right_pw[i,j,k] = pw[i+1,j,k] - pw[i-1,j,k]
rainc = data["rainc"]
alpha = np.zeros((49,180,360))
for i in range(1,47):
    for j in range(0,179):
        for k in range(0,359):
            if right_pw[i,j,k] != 0:
                alpha[i,j,k] = rainc[i,j,k] / pw[i,j,k]
np.save("alpha",alpha)