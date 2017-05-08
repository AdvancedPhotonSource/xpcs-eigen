# -*- coding: utf-8 -*-
# @Author: Faisal Khan
# @Date:   2016-11-15 14:57:56
# @Last Modified by:   Faisal Khan
# @Last Modified time: 2017-01-19 16:37:31

import sys

sys.path.append("io")
sys.path.append("filter")

from imm_file import *
from pixel_sum import *
from frame_sum import *
import cPickle as pickle

imm = IMM("/local/fkhan/data/xpcs/A006_PinarFe5wtpPS_7p35keV_20H20V_USID_Sq1_001_00001-00266.imm")
arr = imm.array()

print(arr.shape)
print (arr.dtype)

print (arr[0:9, 0:20])
print (arr[1:10, 0:20])
#print (np.sum(arr, axis=0))