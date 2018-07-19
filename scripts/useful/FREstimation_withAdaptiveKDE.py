
# coding: utf-8

# Here is code used to explore the use of https://github.com/cooperlab/AdaptiveKDE

# In[28]:


from GJEphys.ssvkernel import ssvkernel
import pandas as pd
import matplotlib.pyplot as plt
from GJEphys.matplotlibRCParams import mplPars
import seaborn as sns
import os
from GJEphys.pdColumnNameMapCont import mdFN
import numpy as np


# In[22]:


dataFile = os.path.expanduser("~/DataAndResults/ephys/Results/ContStimPSTH/data_15Nrns.xlsx")


# In[36]:


expID = "130313-4Rh"
trialName = "Trial5"


# In[24]:


dataDF = pd.read_excel(dataFile)


# In[37]:


trialCrit = np.logical_and(dataDF[mdFN["expID"]] == expID, dataDF[mdFN["trialName"]] == trialName)
trialData = dataDF.loc[trialCrit, :]
trialSpikeTimes = trialData["Spike Time (s)"]
# print(trialSpikeTimes)


# In[44]:


tEst = np.arange(-0.5, 1.5, 0.02)


# In[43]:


res = ssvkernel(x=trialSpikeTimes.values, tin=tEst)

