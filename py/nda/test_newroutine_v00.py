import srn_nda_newroutine_jup as nda
import glob

file_list = glob.glob('/Users/baptiste/Projets/ObsNancay/NDA/NewRoutine/*.dat')

for item in file_list:
    nda.build_edr_cdf(item, 'config.json', True)
