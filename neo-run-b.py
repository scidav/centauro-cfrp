# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__


from neolab import *
import os
import glob
import numpy as np




# ## PLEASE DO NOT DELETE THE FOLLOWING LINES ==============================================================
# def getValues(keyStr):
#     databases = [v for v in glob.glob("{0:}/*.odb".format(wd)) if keyStr in v.split("\\")[-1]]
#     for v in databases: 
#         extractCurves(v)
#     #
#     getMax()

# ## PLEASE DO NOT DELETE THE UPSIDE LINES ==============================================================
keyStr = "M"

databases = [v for v in glob.glob("{0:}/*.odb".format(wd)) if keyStr in v.split("\\")[-1]]
for v in databases: 
    extractCurves(v)


getMax()


# # Uncomment for get Hashin history (REPAIR models valid only)
# getHashin()
# exportLoadCurves()









