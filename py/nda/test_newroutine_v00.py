import srn_nda_newroutine_jup as nda
import glob

file_list = glob.glob('/Users/baptiste/Projets/ObsNancay/NDA/NewRoutine/data/*.dat')

for item in file_list:
    header,var = nda.load_dat_data(item, 1000)
