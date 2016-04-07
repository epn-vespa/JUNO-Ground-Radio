import numpy as np
import os
import os.path
import datetime

def load_rt1_data(file_rt1):
	
	print 'Loading RT1 file [NDA/RoutineJup]: '+file_rt1
	record_size = 405
	file_size = os.path.getsize(file_rt1)
	nspectra = file_size / record_size - 1
	
	print 'Record size = '+str(record_size)+' Bytes'
	print 'File size   = '+str(file_size)+' Bytes'
	print 'Nb of Spectra = '+str(nspectra)
	
	dt_raw = np.dtype([('ss','a',405),('bb',[('rec_hr','b'),('rec_min','b'),('rec_sec','b'),('rec_ms','b'),('data','b',400),('status','b')],nspectra)])
	raw = np.fromfile(file_rt1, dtype=dt_raw)
	
	print 'Interpreting Header...'
	ss = raw['ss'][0]
	
	if (int(ss[56:58]) < 90): 
		century_str1 = '20'
	else:
		century_str1 = '19'
	if (int(ss[47:49]) < 90):
		century_str2 = '20' 
	else:
		century_str2 = '19'
	
	header = {}
	header['fmin']    = float(ss[1:3]) 
	header['fmax']    = float(ss[3:5])
	header['fstep']   = float(ss[5:8])
	header['reflev']  = float(ss[8:11])
	header['time_res'] = float(ss[11:16])
	header['power_res'] = float(ss[16:18])
	header['rf_filter0'] = int(ss[22])
	header['rf_filter1'] = int(ss[28])
	
	if (header['rf_filter1'] == 0): 
		header['rf_filter_time1'] = datetime.date.min
	else:
		header['rf_filter_time1'] = datetime.datetime(int(century_str1+ss[56:58]),int(ss[53:55]),int(ss[50:52]),int(ss[29:31]),int(ss[32:34]),0)
		if (int(ss[29:31]) <= int(ss[23:25])):
			header['rf_filter_time1'] += datetime.timedelta(days=1)
	
	header['rf_filter2'] = int(ss[34])
	
	if (header['rf_filter2'] == 0):
		header['rf_filter_time2'] = datetime.date.min
	else: 
		header['rf_filter_time2'] = datetime.datetime(int(century_str1+ss[56:58]),int(ss[53:55]),int(ss[50:52]),int(ss[35:37]),int(ss[38:40]),0)
		if (int(ss[35:37]) <= int(ss[23:25])):
			header['rf_filter_time2'] += datetime.timedelta(days=1)
	
	header['meridian_time'] = datetime.datetime(int(century_str1+ss[47:49]),int(ss[44:46]),int(ss[41:43]),int(ss[18:20]),int(ss[20:22]),0)
	header['start_time']    = datetime.datetime(int(century_str1+ss[56:58]),int(ss[53:55]),int(ss[50:52]),int(ss[23:25]),int(ss[26:28]),0)
	header['stop_time']     = datetime.datetime(int(century_str1+ss[56:58]),int(ss[53:55]),int(ss[50:52]),int(ss[59:61]),int(ss[62:64]),0)
	
	print 'Loading data...'
	
	base_date = datetime.datetime(header['start_time'].year,header['start_time'].month,header['start_time'].day)
	data_time = np.array([base_date + datetime.timedelta(hours=int(raw['bb']['rec_hr'][0][i]),minutes=int(raw['bb']['rec_min'][0][i]),seconds=int(raw['bb']['rec_sec'][0][i]),milliseconds=int(raw['bb']['rec_ms'][0][i])*10) for i in np.arange(nspectra)])
	
	for i in np.arange(nspectra):
		if (data_time[i] < data_time[0]):
			data_time[i] += datetime.timedelta(days=1)
	
	data = raw['bb']['data'][0]
	data_status = raw['bb']['status'][0]
	
	return header, data_time.tolist(), data, data_status.tolist()

def build_edr_data(file):
	
	print 'Loading data...'
	header, date, spectra, status_code = load_rt1_data(file)
	
	nspectra = len(date)
	
	print 'Building polarization indices...'
	index_lh = [((ii % 2) == 1) for ii in range(0, nspectra)]
	index_rh = [ (not ii) for ii in index_lh]
	cnt_lh = sum(index_lh)
	cnt_rh = sum(index_rh)
	
	print 'LH data :'+str(cnt_lh)+' sweeps'
	print 'RH data :'+str(cnt_rh)+' sweeps'
	
	out = {}
	out['header'] = header
	out['datetime'] = date
	out['spectra'] = spectra
	out['index_lh'] = index_lh
	out['index_rh'] = index_rh
	out['status'] = status_code
	return out

def load_strpou_data(file_pou):
	
	print 'Loading pointing file: '+file_pou
	f = open(file_pou, 'r')
	raw = f.read()
	f.close()
	print 'Splitting raw data into lines'
	lines = raw.splitlines()
	print 'Defining PointingTable array' 
	PointingTable = []
	i = 0
	for line in lines:
		if (line[0] != '#'):
			PointingTable.append({})
			if (int(line[6:8]) < 90):
				century = 2000
			else:
				century = 1900
			PointingTable[i]['datetime'] = datetime.datetime(century+int(line[6:8]),int(line[3:5]),int(line[0:2]),int(line[9:11]),int(line[12:14]),int(line[15:17]))
			PointingTable[i]['phase_RA'] = []
			PointingTable[i]['phase_LA'] = []
			for j in range(0,8):
				PointingTable[i]['phase_RA'].append(int(line[18+j]))
				PointingTable[i]['phase_LA'].append(int(line[48+j]))
			PointingTable[i]['delay_EW_RA'] = int(line[27:30])
			PointingTable[i]['delay_EW_LA'] = int(line[57:60])
			PointingTable[i]['delay_NS_RA'] = int(line[31:34])
			PointingTable[i]['delay_NS_LA'] = int(line[61:64])
			PointingTable[i]['Field_filter_RA'] = int(line[35])
			PointingTable[i]['Field_filter_LA'] = int(line[55])
			PointingTable[i]['Lab_filter_RA'] = []
			PointingTable[i]['Lab_filter_LA'] = []
			for j in range(0,4):
				PointingTable[i]['Lab_filter_RA'].append(int(line[37+j*2]))
				PointingTable[i]['Lab_filter_LA'].append(int(line[67+j*2]))
			if (line[45:47] == 'OF'):
				PointingTable[i]['noise_at_RA'] = 255
			else:
				PointingTable[i]['noise_at_RA'] = int(line[45:47])
			if (line[75:77] == 'OF'):
				PointingTable[i]['noise_at_LA'] = 255
			else:
				PointingTable[i]['noise_at_LA'] = int(line[75:77])
			PointingTable[i]['Interarray_delay'] = int(line[78:81])
			i += 1
	
	cnt_pou = i
	return PointingTable

def rt1_to_ddr_calib(file_data,file_pointing,ignore_first_cal):
	
	edr = build_edr_data(file_data)
	out = edr
	#header, date, data, status_code, index_lh, index_rh
	pointing = load_strpou_data(file_pointing)
	
	cnt_lh = sum(edr['index_lh'])
	cnt_rh = sum(edr['index_rh'])
	
	print 'LH data :'+str(cnt_lh)+' sweeps'
	print 'RH data :'+str(cnt_rh)+' sweeps'
	
	print 'Calibration sequence identification...'
	index_cal = [(x['noise_at_LA'] != 255) | (x['noise_at_RA']!= 255) for x in pointing] 
	cnt_cal = sum([int(i) for i in index_cal])
	print 'Found '+str(cnt_cal)+' calibration intervals.'
	
	cal_sequences = {}
	cal_sequences['start_time'] = []
	cal_sequences['stop_time'] = []
	for x in pointing:
		if x['noise_at_LA'] == 30:
			cal_sequences['start_time'].append({})
			cal_sequences['start_time']['']
		cal_sequences[int(i)]
	
	if (ignore_first_cal == True):
		print "IGNORE_FIRST_CAL: not yet implemented"
	
	return out


	
	

