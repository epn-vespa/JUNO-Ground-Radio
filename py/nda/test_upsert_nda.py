#! /usr/bin/python
import os
import glob
from spacepy import pycdf

data_dir = "/databf/nda/web/jupiter/data/"
plot_dir = "/databf/nda/web/jupiter/ql/B-W/"
data_yr = "2017"
data_mo = "01"
cdf_list = sorted(glob.glob(os.path.join(data_dir, data_yr, data_mo, "*.cdf")))
pdf_list = sorted(glob.glob(os.path.join(plot_dir, data_yr, data_mo, "J??????.pdf")))
ncdf = len(cdf_list)
npdf = len(pdf_list)

if ncdf != npdf:
    print 'error !'
else:
    for i in range(ncdf):
        os.system('python upsert_nda.py --verb {} {}'.format(cdf_list[i], pdf_list[i]))