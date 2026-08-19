import numpy
import rpy2
import pandas as pd
import os
from vectra_lib_v3 import nearest_neighbor

output=extract_data('/Users/ee699/working/TRIC/HIPC/Human_LN-Healthy/whole-panel_lung_ln/out-phenomenalist/D437-IRF8-segmentation/sub-clustering/pcf-inputs/')
nn=nearest_neighbor(output)
