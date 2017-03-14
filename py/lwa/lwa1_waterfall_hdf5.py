import numpy as np
import datetime
import os
import os.path
from astropy.time import Time as astroTime
from spacepy import pycdf
import h5py

def get_paths():
	master_path = './'
	cdfbin_path = '/Applications/cdf/cdf36_0-dist/bin/'
	cdfout_path = './'
	return 	master_path, cdfbin_path, cdfout_path

def master_cdf_name(obsty_id, instr_id, level, cdf_version):
	return '{}_{}_{}_000000000000_000000000000_V{}.cdf'.format(obsty_id, instr_id, level, cdf_version)

def master_skt_name(obsty_id, instr_id, level, cdf_version):
	return '{}_{}_{}_000000000000_000000000000_V{}.skt'.format(obsty_id, instr_id, level, cdf_version)

class data:

	def __init__(self,file):

		# test file = '055917_000007389_Jupiter-waterfall-low.hdf5'
		f_hdf = h5py.File(file,'r')
	
		ndata = len(f_hdf['Observation1']['time'][:])
		nfreq = len(f_hdf['Observation1']['Tuning1']['freq'][:])
	
		edr = {}
	
		edr['header'] = {}
		edr['header']['obsty_id'] = 'lwa1'
		edr['header']['instr_id'] = 'waterfall_low'
		edr['header']['frequency'] = f_hdf['Observation1']['Tuning1']['freq'][:]
			
		edr['datetime'] = []	
		base_date = datetime.datetime(1970, 1, 1)
		for t in f_hdf['Observation1']['time'][...]:
			edr['datetime'].append(base_date + datetime.timedelta(seconds=t))
	
		edr['flux_XX'] = f_hdf['Observation1']['Tuning1']['XX'][:]
		edr['flux_YY'] = f_hdf['Observation1']['Tuning1']['YY'][:]
		edr['saturation'] = f_hdf['Observation1']['Tuning1']['Saturation'][:]
	
		self.edr=edr

	def to_cdf(self):

		edr = self.edr 
		
		dat_version = '01'
		sft_version = '01'
		cdf_version = '01'
		
		master_path, cdfbin_path, cdfout_path = get_paths()
	
		master_cdf = master_cdf_name(edr['header']['obsty_id'], edr['header']['instr_id'], 'edr', cdf_version)
		skelet_cdf = master_skt_name(edr['header']['obsty_id'], edr['header']['instr_id'], 'edr', cdf_version)
	
	
		if os.path.exists(master_path+master_cdf):
			os.remove(master_path+master_cdf)
		os.system(cdfbin_path+'skeletoncdf -cdf '+master_path+master_cdf+' '+master_path+skelet_cdf)
		
		ndata = len(edr['datetime'])
		nfreq = len(edr['header']['frequency'])
		
		jul_date = astroTime(edr['datetime'],format="datetime",scale="utc").jd.tolist()
		
		frequency = edr['header']['frequency'] 
	
		f_cdf_name = "{}_{}_edr_{:%Y%m%d%H%M}_{:%Y%m%d%H%M}_V{}.cdf".format(edr['header']['obsty_id'], edr['header']['instr_id'], edr['datetime'][0], edr['datetime'][ndata-1],cdf_version)
		if os.path.exists(cdfout_path+f_cdf_name):
			os.remove(cdfout_path+f_cdf_name)
		f_cdf = pycdf.CDF(cdfout_path+f_cdf_name, master_path+master_cdf)
	
		# SETTING PDS GLOBAL ATTRIBUTES
		f_cdf.attrs['PDS_Observation_start_time'] = edr['datetime'][0].isoformat()+'Z'
		f_cdf.attrs['PDS_Observation_stop_time'] = edr['datetime'][ndata-1].isoformat()+'Z'
	
		# SETTING VESPA GLOBAL ATTRIBUTES
		f_cdf.attrs['VESPA_time_min'] = jul_date[0]
		f_cdf.attrs['VESPA_time_max'] = jul_date[ndata-1]
		f_cdf.attrs['VESPA_time_sampling_step_min'] = np.amin([jul_date[i+1]-jul_date[i] for i in range(0,ndata-2)])*86400.
		f_cdf.attrs['VESPA_time_sampling_step_max'] = np.amax([jul_date[i+1]-jul_date[i] for i in range(0,ndata-2)])*86400.
	
		f_cdf.attrs['VESPA_spectral_range_min']  = np.amin(frequency)
		f_cdf.attrs['VESPA_spectral_range_max']  = np.amax(frequency)
		f_cdf.attrs['VESPA_spectral_sampling_step_min'] = np.amin([frequency[i+1]-frequency[i] for i in range(0,nfreq-2)])
		f_cdf.attrs['VESPA_spectral_sampling_step_max'] = np.amax([frequency[i+1]-frequency[i] for i in range(0,nfreq-2)])
		f_cdf.attrs['VESPA_spectral_resolution_min'] = f_cdf.attrs['VESPA_spectral_sampling_step_min']
		f_cdf.attrs['VESPA_spectral_resolution_max'] = f_cdf.attrs['VESPA_spectral_sampling_step_max']
	
		# SETTING OTHER GLOBAL ATTRIBUTES
		f_cdf.attrs['Logical_file_id'] = f_cdf_name
		f_cdf.attrs['Data_version'] = dat_version
		f_cdf.attrs['Skeleton_version'] = cdf_version
		f_cdf.attrs['Software_version'] = sft_version
		f_cdf.attrs['Software_language'] = 'python'
	
		# SETTING VARIABLES
		f_cdf['Epoch'] = edr['datetime']
		f_cdf['ISO_DATE'] = [d.isoformat()+'Z' for d in edr['datetime']]
		f_cdf['JD_TIME'] = jul_date
		f_cdf['FLUX_XX'] = edr['flux_XX']
		f_cdf['FLUX_YY'] = edr['flux_YY']
		f_cdf['SATURATION'] = edr['saturation']
		f_cdf['Frequency'] = frequency
	
		date_start = edr['datetime'][0]
		date_stop = edr['datetime'][ndata-1]
		date_start_round = edr['datetime'][0].replace(minute=0, second=0, microsecond=0)
		date_stop_round = edr['datetime'][ndata-1].replace(minute=0, second=0, microsecond=0)+datetime.timedelta(hours=1)
	
		# SETTING VARIABLES ATTRIBUTES
		f_cdf['Epoch'].attrs['VALIDMIN'] = date_start
		f_cdf['Epoch'].attrs['VALIDMAX'] = date_stop
		f_cdf['Epoch'].attrs['SCALEMIN'] = date_start_round
		f_cdf['Epoch'].attrs['SCALEMAX'] = date_stop_round
	
		f_cdf['ISO_DATE'].attrs['VALIDMIN'] = date_start.isoformat()+'Z'
		f_cdf['ISO_DATE'].attrs['VALIDMAX'] = date_stop.isoformat()+'Z'
		f_cdf['ISO_DATE'].attrs['SCALEMIN'] = date_start_round.isoformat()+'Z'
		f_cdf['ISO_DATE'].attrs['SCALEMAX'] = date_stop_round.isoformat()+'Z'
	
		f_cdf['JD_TIME'].attrs['VALIDMIN'] = astroTime(date_start,format="datetime",scale="utc").jd
		f_cdf['JD_TIME'].attrs['VALIDMAX'] = astroTime(date_stop,format="datetime",scale="utc").jd
		f_cdf['JD_TIME'].attrs['SCALEMIN'] = astroTime(date_start_round,format="datetime",scale="utc").jd
		f_cdf['JD_TIME'].attrs['SCALEMAX'] = astroTime(date_stop_round,format="datetime",scale="utc").jd
	
		f_cdf.close()
	
