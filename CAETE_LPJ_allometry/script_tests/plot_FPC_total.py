import pandas as pd 

import matplotlib.pyplot as plt 

FPC = pd.read_csv('totalFPC.csv')

carbon = pd.read_csv('carbons.csv')
carbon.columns = ['cleaf','cwood','croot','pls','time']
leaf = carbon.cleaf
wood = carbon.cwood
root = carbon.croot
pls = carbon.pls
time = carbon.time

# plt.plot(time,leaf)

# plt.show()

plt.plot(FPC)

plt.show()