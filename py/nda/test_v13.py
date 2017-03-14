import srn_nda_routine_jup as nda
import glob

file_list = glob.glob('/Users/baptiste/Projets/ObsNancay/NDA/RT1/headers/data/*/*.RT1')

for file in file_list: 
	nda.build_edr_cdf(file,True)


