import numpy as np
import os
import os.path
import struct
import datetime
from spacepy import pycdf
import json


def get_paths(config_file, cdf_version):
    """
    Sets up the paths for processing.
    :param cdf_version: version of CDF dataset
    :param config_file: configuration file (json formatted)
    :return master_path: Path to master CDF
    :return cdfbin_path: Path to GSFC CDF library
    :return cdfout_path: Path to output directory
    """

    with open(config_file) as fh:
        config = json.load(fh)

    return config['path']['master_path'], config['path']['cdfbin_path'], config['path']['cdfout_path']


def get_versions(config_file):
    """
    Sets up the versions for processing.
    :param config_file: configuration file (json formatted)
    :return cdf_version: Version of CDF skeleton
    :return sft_version: Version of current software
    :return dat_version: Version of input data
    """

    with open(config_file) as fh:
        config = json.load(fh)

    return config['version']['cdf'], config['version']['sft'], config['version']['dat']


def master_cdf_name(obsty_id, instr_id, recvr_id, level, cdf_version):
    """
    Builds the CDF master file name
    :param obsty_id: Observatory ID (short name)
    :param instr_id: Instrument ID (short name)
    :param recvr_id: Receiver ID (short name)
    :param level: data level (edr or ddr)
    :param cdf_version: dataset version
    :return : name of the CDF master
    """

    return '{}_{}_{}_{}_000000000000_000000000000_V{}.cdf'.format(obsty_id, instr_id, recvr_id, level, cdf_version)


def master_skt_name(obsty_id, instr_id, recvr_id, level, cdf_version):
    """
    Builds the CDF skeleton file name
    :param obsty_id: Observatory ID (short name)
    :param instr_id: Instrument ID (short name)
    :param recvr_id: Receiver ID (short name)
    :param level: data level (edr or ddr)
    :param cdf_version: dataset version
    :return : name of the CDF skeleton
    """

    return '{}_{}_{}_{}_000000000000_000000000000_V{}.skt'.format(obsty_id, instr_id, recvr_id, level, cdf_version)


def decode_header(hdr_raw):
    """
    Decodes the DAT file header (see SRN NDA NEW ROUTINE JUP documentation)
    :param hdr_raw: Raw header data
    :return header: Header data dictionary
    """

    hdr_fmt = '<68I1l2048f2048l'
    hdr_val = struct.unpack(hdr_fmt, hdr_raw)

    header = dict()
    header['siz_tet'] = hdr_val[0]       # Header Size
    header['sel_prod0'] = hdr_val[1]     # selected products 0
    header['sel_prod1'] = hdr_val[2]     # selected products 1
    header['acc'] = hdr_val[3]           # accumulating factor
    header['subband'] = hdr_val[4:68]    # selected sub-bands
    header['nfft'] = hdr_val[68]         # Number of FFT points
    header['frq'] = hdr_val[69:2117]     # frequency values
    header['ifrq'] = hdr_val[2117:4165]  # frequency indices

    return header


class JDateTime:

    def __init__(self, jdatetime):
        self.jdatetime = jdatetime

    def to_jd(self):
        return float(self.jdatetime[0]) + \
               (float(self.jdatetime[1]) + float(self.jdatetime[2]) / float(self.jdatetime[3])) / 86400.

    def to_datetime(self):
        return datetime.datetime(2000, 1, 1, 12)-datetime.timedelta(2451545-self.to_jd())


def load_dat_data(input_file):
    """
    Loads RT1 data
    :param input_file: input RT1 file
    :return:
    """

    print '** Loading DAT file [SRN/NDA/NewRoutine]: '+input_file
    header_size = 16660
    record_size = 32832
    file_size = os.path.getsize(input_file)
    nspectra = (file_size - header_size)/record_size

    print '   Header size = '+str(header_size)+' Bytes'
    print '   Record size = '+str(record_size)+' Bytes'
    print '   File size   = '+str(file_size)+' Bytes'
    print '   Nb of Spectra = '+str(nspectra)

    # There are various .DAT file naming conventions (but we don't really care :-)
    # - from Feb 2017: JYYYYMMDD_HHMMSS_YYYYMMDD_HHMMSS_Rou.dat,
    #   with start and stop times (stored /datadf/nda/web/jupiter/data/YYYY/MM)
    # - during Jan 2017: JYYYYMMDD_HHMMSS_Rou.dat, with start time
    #   (stored /datadf/nda/web/jupiter/data/YYYY/MM)
    # - from 2012 to Dec 2015: JYYYYMMDD_HHMMSS_YYYYMMDD_HHMMSS.dat,
    #   with start and stop times (stored /databf/nda/NewRoutine/YYYY/)
    # - during 2016: JYYYYMMDD_HHMMSS_YYYYMMDD_HHMMSS.dat or JYYYYMMDD_HHMMSS.dat
    #   with start and stop times (stored /databf/nda/NewRoutine/2016/)

    ff = open(input_file, 'rb')
    print '** Loading Header from file.'
    header = decode_header(ff.read(header_size))

    sel_chan = [int(ii) for ii in
                list('{:032b}'.format(header['sel_prod0'])[::-1] + '{:032b}'.format(header['sel_prod1'])[::-1])]
    nbchan = sum(sel_chan)

    print '** Loading Data from file.'

    record_fmt = '<8I'
    for iichan in range(0, nbchan):
        record_fmt = '{}{}'.format(record_fmt, '2I2048f')

    packet_size = 10000

    data = list()
    jdatetime = list()

    for ii in range(0, nspectra, packet_size):
        read_size = packet_size
        if (ii+packet_size)*record_size+header_size > file_size:
            read_size = nspectra-ii
        print '   {:05d} to {:05d}: reading {} records'.format(ii,ii+read_size,read_size)
        packet_raw = ff.read(record_size*read_size)
        for jj in range(0, read_size):
            rec_val = struct.unpack(record_fmt, packet_raw[jj*record_size:(jj+1)*record_size])
            record = dict()
            record['header'] = dict()
            record['header']['magic'] = rec_val[0]
            record['header']['cnt'] = rec_val[1]
            record['header']['jdatetime'] = rec_val[2:6]
            record['data'] = list()
            for ichan in range(0, nbchan):
                record['data'].append(dict())
                record['data'][ichan]['magic'] = rec_val[8+ichan*2050]
                record['data'][ichan]['spec_typ'] = rec_val[9+ichan*2050]
                record['data'][ichan]['spectrum'] = rec_val[10+ichan*2050:2060+ichan*2050]
            data.append(record)
            jdt = JDateTime(record['header']['jdatetime'])
            jdatetime.append(jdt)

    # ichan = 0 ==> Left-Handed Array (LL)
    # ichan = 1 ==> Real part of cross-correlation (CR) [empty before 2013/10/18]
    # ichan = 2 ==> Imaginary part of cross-correlation (CI) [empty before 2013/10/18]
    # ichan = 3 ==> Right-Handed Array (RR)

    # A FAIRE
    # - recalculer l'heure + jour julien
    # - filtrer sur 10 MHZ a 40 MHz (verifier filtre labos ? )

    var = dict()
    var['data'] = data
    var['jdatetime'] = jdatetime

    header['obsty_id'] = 'srn'
    header['instr_id'] = 'nda'
    header['recvr_id'] = 'newroutine'

    # header['freq_min'] = float(raw_header['freq_min']) 	# MHz
    # header['freq_max'] = float(raw_header['freq_max']) 	# MHz
    # header['freq_res'] = float(raw_header['freq_res'])  # kHz
    # header['frequncy'] = np.arange(400)/400.*(header['freq_max']-header['freq_min'])+header['freq_min']
    # header['freq_stp'] = (header['freq_max']-header['freq_min'])/0.4   # kHz
    # header['time_stp'] = 1.
    #
    # header['ref_levl'] = float(raw_header['ref_levl'])  # dBm
    # header['swp_time'] = float(raw_header['swp_time'])  # ms
    # header['time_exp'] = 0.875        					# ms
    # header['powr_res'] = float(raw_header['powr_res'])
    #
    # header['raw_header_version'] = header_version
    #
    # # Meridian Date-Time crossing
    # if header_version == 1:
    #     # NB: For header version 1, the HH:MM of the meridian crossing is not provided.
    #     header['meridian_datetime'] = datetime.datetime(int(raw_header['merid_yr']), int(raw_header['merid_mo']),
    #                                                     int(raw_header['merid_dd']), 0, 0, 0)
    # else:
    #     header['meridian_datetime'] = datetime.datetime(int(raw_header['merid_yr']), int(raw_header['merid_mo']),
    #                                                     int(raw_header['merid_dd']), int(raw_header['merid_hh']),
    #                                                     int(raw_header['merid_mm']), 0)
    #
    # # Extracting observation Start and Stop Times from data
    # raw_header['beg_hh'] = int(bb['rec_hr'][0])
    # raw_header['beg_mm'] = int(bb['rec_min'][0])
    # raw_header['beg_ss'] = int(bb['rec_sec'][0])
    # raw_header['beg_ms'] = int(bb['rec_cs'][0])*10
    #
    # raw_header['end_hh'] = int(bb['rec_hr'][-1])
    # raw_header['end_mm'] = int(bb['rec_min'][-1])
    # raw_header['end_ss'] = int(bb['rec_sec'][-1])
    # raw_header['end_ms'] = int(bb['rec_cs'][-1])*10
    #
    # # Building Observation Start and Stop Date-Times
    # # Assuming observation dates same as meridian dates (it is fixed below if required)
    # header['start_datetime'] = \
    #     datetime.datetime(int(raw_header['merid_yr']), int(raw_header['merid_mo']), int(raw_header['merid_dd']),
    #                       raw_header['beg_hh'], raw_header['beg_mm'], raw_header['beg_ss'],
    #                       raw_header['beg_ms']*1000)
    # header['stop_datetime'] = \
    #     datetime.datetime(int(raw_header['merid_yr']), int(raw_header['merid_mo']), int(raw_header['merid_dd']),
    #                       raw_header['end_hh'], raw_header['end_mm'], raw_header['end_ss'],
    #                       raw_header['end_ms']*1000) + \
    #     datetime.timedelta(milliseconds=header['swp_time'])
    # header['start_time'] = header['start_datetime'].time()
    # header['stop__time'] = header['stop_datetime'].time()
    # # Fixing start and stop dates when midnight occurs during observation
    # if header['start_time'] > header['stop__time']:
    #     # we don't check header_version = 1 as the corresponding data do not fall into this case (pfiou...)
    #     if header['meridian_datetime'].time() > header['start_time']:
    #         # if meridian time is between start time and midnight, then add 1 day to stop time
    #         header['stop_datetime'] += datetime.timedelta(days=1)
    #     else:
    #         # if meridian time is between midnight and stop time, then subtract 1 day to start time
    #         header['start_datetime'] += datetime.timedelta(days=-1)
    #
    # # Checking header Observation Stop Time, when it exist (header version 6)
    # if header_version == 6:
    #     if int(raw_header['h_stp_hh']) != header['stop__time'].hour:
    #         print 'Warning. Header Stop time (HH) inconsistent with data records!'
    #         print '         Keeping data record value. Header = {:02d}:{:02d} / Data = {:02d}:{:02d}'\
    #             .format(int(raw_header['h_stp_hh']), int(raw_header['h_stp_mm']),
    #                     header['stop__time'].hour, header['stop__time'].minute)
    #     if int(raw_header['h_stp_mm']) != header['stop__time'].minute:
    #         print 'Warning. Header Stop time (MM) inconsistent with data records!'
    #         print '         Keeping data record value. Header = {:02d}:{:02d} / Data = {:02d}:{:02d}'\
    #             .format(int(raw_header['h_stp_hh']), int(raw_header['h_stp_mm']),
    #                     header['stop__time'].hour, header['stop__time'].minute)
    #
    # # RF Filter configuration
    # if header_version < 5:
    #     # NB: For header version 1 to 4, the RF Filter configuration is not provided (see .POU or .COM files)
    #     header['rf_filt0'] = 0
    #     header['rf_filt1'] = 0
    #     header['rf_filt2'] = 0
    #
    #     header['rf_filt0_time'] = datetime.date.min
    #     header['rf_filt1_time'] = datetime.date.min
    #     header['rf_filt2_time'] = datetime.date.min
    # else:
    #     header['rf_filt0'] = int(raw_header['rf0_sele'])
    #     header['rf_filt1'] = int(raw_header['rf1_sele'])
    #     header['rf_filt2'] = int(raw_header['rf2_sele'])
    #
    #     header['rf_filt0_time'] = header['start_datetime']
    #     if header['rf_filt1'] == 0:
    #         header['rf_filt1_time'] = datetime.date.min
    #     else:
    #         header['rf_filt1_time'] = \
    #             datetime.datetime.combine(header['start_datetime'].date(),
    #                                       datetime.time(int(raw_header['rf1_hour']), int(raw_header['rf1_minu'])))
    #         if int(raw_header['rf1_hour']) < int(raw_header['beg_hh']):
    #             header['rf_filt1_time'] += datetime.timedelta(days=1)
    #
    #     if header['rf_filt2'] == 0:
    #         header['rf_filt2_time'] = datetime.date.min
    #     else:
    #         header['rf_filt2_time'] = \
    #             datetime.datetime.combine(header['start_datetime'].date(),
    #                                       datetime.time(int(raw_header['rf2_hour']), int(raw_header['rf2_minu'])))
    #         if int(raw_header['rf2_hour']) < int(raw_header['beg_hh']):
    #             header['rf_filt2_time'] += datetime.timedelta(days=1)
    #
    # print 'Loading data...'
    #
    # var = dict()
    # var['time'] = np.array(
    #     [datetime.datetime.combine(header['start_datetime'].date(),
    #                                datetime.time(int(bb['rec_hr'][i]), int(bb['rec_min'][i]),
    #                                              int(bb['rec_sec'][i]), int(bb['rec_cs'][i])*10000))
    #      for i in np.arange(nspectra)]).tolist()
    #
    # for i in np.arange(nspectra):
    #     if var['time'][i] < var['time'][0]:
    #         var['time'][i] += datetime.timedelta(days=1)
    #
    # var['data'] = bb['data']
    # var['status'] = bb['status'].tolist()
    # var['offset'] = np.array([float("{0:.2f}".format(0.875*i)) for i in np.arange(400)])

    return header, var


# def build_edr_data(input_file):
#     """
#     Builds EDR data file from RT2 file
#     :param input_file: input RT1 file
#     :return:
#     """
#
#     print 'Loading data...'
#     header, var = load_rt1_data(input_file)
#
#     nspectra = len(var['time'])
#
#     print 'Building polarization indices...'
#     index_lh = [((ii % 2) == 1) for ii in range(0, nspectra)]
#     index_rh = [(not ii) for ii in index_lh]
#     cnt_lh = sum(index_lh)
#     cnt_rh = sum(index_rh)
#
#     print 'LH data :'+str(cnt_lh)+' sweeps'
#     print 'RH data :'+str(cnt_rh)+' sweeps'
#
#     out = dict()
#     out['header'] = header
#     out['time'] = list()
#     out['ll_flux'] = list()
#     out['rr_flux'] = list()
#     out['status'] = list()
#     out['rr_sweep_time_offset'] = list()
#
#     for ii in np.arange(nspectra/2):
#         out['time'].append(var['time'][ii*2])
#         out['ll_flux'].append(var['data'][:][ii*2])
#         out['rr_flux'].append(var['data'][:][ii*2+1])
#         tmp_status = list()
#         tmp_status.append(var['status'][ii*2])
#         tmp_status.append(var['status'][ii*2+1])
#         out['status'].append(tmp_status)
#         rr_sweep_time_offset = var['time'][ii*2+1]-var['time'][ii*2]
#         out['rr_sweep_time_offset'].append(rr_sweep_time_offset.total_seconds())
#
#     out['sweep_time_offset_ramp'] = var['offset']
#
# #   Calibration sweep detection [broken, so removed from processing for the current version]
# #   out['cal_start_time'] = []
# #   out['cal_stop_time'] = []
# #   out['cal_attenuation'] = []
# #   cal_mode = -1
# #   for ii in np.arange(nspectra):
# #       if var['status'][ii] == 17:
# #           out['cal_start_time'].append(var['time'][ii+1])
# #           if cal_mode == -1:
# #               out['cal_attenuation'].append(30.)
# #               cal_mode = 0
# # 	        elif cal_mode == 0:
# # 	            out['cal_stop_time'].append(var['time'][ii])
# #               out['cal_attenuation'].append(20.)
# #               cal_mode = 1
# # 			elif cal_mode == 1:
# # 				out['cal_stop_time'].append(var['time'][ii])
# # 				out['cal_attenuation'].append(10.)
# # 				cal_mode = 2
# # 			elif cal_mode == 2:
# # 				out['cal_stop_time'].append(var['time'][ii])
# # 				out['cal_attenuation'].append(00.)
# # 				cal_mode = -2
# # 		if var['status'][ii] == 0:
# # 			if cal_mode == -2:
# # 				out['cal_stop_time'].append(var['time'][ii])
# # 				cal_mode = -1
#
#     return out
#
#
# def build_edr_cdf(input_file, config_file='config.json', build_cdf_master=False):
#     """
#     Builds the EDR level CDF file
#     :param input_file: input RT1 file
#     :param config_file: configuration file (default = config.json in same directory)
#     :param build_cdf_master: option to (re)build CDF master (default = False)
#     :return:
#     """
#
#     cdf_version, sft_version, dat_version = get_versions(config_file)
#     master_path, cdfbin_path, cdfout_path = get_paths(config_file, cdf_version)
#
#     edr = build_edr_data(input_file)
#
#     master_cdf = master_cdf_name(edr['header']['obsty_id'], edr['header']['instr_id'], edr['header']['recvr_id'],
#                                  'edr', cdf_version)
#     skelet_cdf = master_skt_name(edr['header']['obsty_id'], edr['header']['instr_id'], edr['header']['recvr_id'],
#                                  'edr', cdf_version)
#
#     master_cdf_path = os.path.join(master_path, master_cdf)
#     skelet_cdf_path = os.path.join(master_path, skelet_cdf)
#
#     if build_cdf_master:
#         if os.path.exists(master_cdf_path):
#             os.remove(master_cdf_path)
#         os.system('{} -cdf {} {}'.format(os.path.join(cdfbin_path,'skeletoncdf'),
#                                          master_cdf_path, skelet_cdf_path))
#
#     print "Raw data file name:"
#     print input_file
#     print "Master CDF file name:"
#     print master_cdf
#
#     ndata = len(edr['time'])
#     frequency = edr['header']['frequncy']
#
#     cdfout_file = "{}_{}_{}_edr_{:%Y%m%d%H%M}_{:%Y%m%d%H%M}_V{}.cdf"\
#         .format(edr['header']['obsty_id'], edr['header']['instr_id'], edr['header']['recvr_id'],
#                 edr['time'][0], edr['time'][ndata-1], cdf_version)
#     cdfout_file_path = os.path.join(cdfout_path, cdfout_file)
#
#     if os.path.exists(cdfout_file_path):
#         os.remove(cdfout_file_path)
#
#     print "Output CDF file name:"
#     print cdfout_file
#
#     cdfout = pycdf.CDF(cdfout_file_path, master_cdf_path)
#
#     # SETTING PDS GLOBAL ATTRIBUTES
#     cdfout.attrs['PDS_Observation_start_time'] = edr['time'][0].isoformat()+'Z'
#     cdfout.attrs['PDS_Observation_stop_time'] = edr['time'][ndata-1].isoformat()+'Z'
#
#     # SETTING VESPA GLOBAL ATTRIBUTES
#     cdfout.attrs['VESPA_time_sampling_step'] = edr['header']['time_stp']        # in seconds
#     cdfout.attrs['VESPA_time_exp'] = edr['header']['time_exp']/1.e3			    # in seconds
#
#     cdfout.attrs['VESPA_spectral_range_min'] = edr['header']['freq_min']*1.e6   # In Hz
#     cdfout.attrs['VESPA_spectral_range_max'] = edr['header']['freq_max']*1.e6   # In Hz
#     cdfout.attrs['VESPA_spectral_sampling_step'] = 75.e3                        # In Hz
#     cdfout.attrs['VESPA_spectral_resolution'] = 30.e3                           # In Hz
#
#     # SETTING SRN-NDA GLOBAL ATTRIBUTES
#
# # 	Calibration sweep detection [broken, so removed from processing for current version]
# # 	tmp_cal_start = []
# # 	tmp_cal_stop  = []
# # 	tmp_cal_atten = []
# # 	for ii in range(len(edr['cal_start_time'])):
# # 		tmp_cal_start.append(edr['cal_start_time'][ii].isoformat()+'Z')
# # 		tmp_cal_stop.append(edr['cal_stop_time'][ii].isoformat()+'Z')
# # 		tmp_cal_atten.append(edr['cal_attenuation'][ii])
# # 	cdfout.attrs['NDA_calibration_start_time'] = tmp_cal_start
# # 	cdfout.attrs['NDA_calibration_stop_time'] = tmp_cal_stop
# # 	cdfout.attrs['NDA_calibration_attenuation'] = tmp_cal_atten
#
#     tmp_rf_filt = []
#     tmp_rf_filt_change = []
#     tmp_rf_filt.append(edr['header']['rf_filt0'])
#     tmp_rf_filt_change.append(edr['header']['start_time'].isoformat()+'Z')
#     if edr['header']['rf_filt1'] != 0:
#         tmp_rf_filt.append(edr['header']['rf_filt1'])
#         tmp_rf_filt_change.append(edr['header']['rf_filt1_time'].isoformat()+'Z')
#     if edr['header']['rf_filt2'] != 0:
#         tmp_rf_filt.append(edr['header']['rf_filt2'])
#         tmp_rf_filt_change.append(edr['header']['rf_filt2_time'].isoformat()+'Z')
#     cdfout.attrs['NDA_rf_filter_selected'] = tmp_rf_filt
#     cdfout.attrs['NDA_rf_filter_time_change'] = tmp_rf_filt_change
#
#     cdfout.attrs['NDA_power_resolution'] = edr['header']['powr_res']
#     cdfout.attrs['NDA_meridian_time'] = edr['header']['meridian_datetime'].isoformat()+'Z'
#     cdfout.attrs['NDA_reference_level'] = edr['header']['ref_levl']
#     cdfout.attrs['NDA_sweep_duration'] = edr['header']['swp_time']/1000.        # in seconds
#     cdfout.attrs['NDA_header_version'] = edr['header']['raw_header_version']
#
#     # SETTING OTHER GLOBAL ATTRIBUTES
#     cdfout.attrs['Logical_file_id'] = cdfout_file
#     cdfout.attrs['Data_version'] = dat_version
#     cdfout.attrs['Skeleton_version'] = cdf_version
#     cdfout.attrs['Software_version'] = sft_version
#     cdfout.attrs['Software_language'] = 'python'
#     cdfout.attrs['Parents'] = os.path.basename(input_file)
#
#     # SETTING VARIABLES
#     cdfout['Epoch'] = edr['time']
#     cdfout['FLUX_RR'] = edr['rr_flux']
#     cdfout['FLUX_LL'] = edr['ll_flux']
#     cdfout['STATUS'] = edr['status']
#     cdfout['SWEEP_TIME_OFFSET_RAMP'] = edr['sweep_time_offset_ramp']
#     cdfout['RR_SWEEP_TIME_OFFSET'] = edr['rr_sweep_time_offset']
#     cdfout['Frequency'] = frequency
#
#     date_start = edr['time'][0]
#     date_stop = edr['time'][ndata-1]
#     date_start_round = edr['time'][0].replace(minute=0, second=0, microsecond=0)
#     date_stop_round = edr['time'][ndata-1].replace(minute=0, second=0, microsecond=0)+datetime.timedelta(hours=1)
#
#     # SETTING VARIABLES ATTRIBUTES
#     cdfout['Epoch'].attrs['VALIDMIN'] = date_start
#     cdfout['Epoch'].attrs['VALIDMAX'] = date_stop
#     cdfout['Epoch'].attrs['SCALEMIN'] = date_start_round
#     cdfout['Epoch'].attrs['SCALEMAX'] = date_stop_round
#
#     cdfout['Frequency'].attrs['SCALEMIN'] = edr['header']['freq_min']
#     cdfout['Frequency'].attrs['SCALEMAX'] = edr['header']['freq_max']
#     cdfout['Frequency'].attrs['UNITS'] = 'MHz'
#
#     cdfout.close()
#
#     print "CDF file written:"
#     print cdfout_file_path


# def load_strpou_data(file_pou):
#
#     print 'Loading pointing file: '+file_pou
#     f = open(file_pou, 'r')
#     raw = f.read()
#     f.close()
#     print 'Splitting raw data into lines'
#     lines = raw.splitlines()
#     print 'Defining pointing_table array'
#     pointing_table = []
#     i = 0
#     for line in lines:
#         if line[0] != '#':
#             pointing_table.append({})
#             if int(line[6:8]) < 90:
#                 century = 2000
#             else:
#                 century = 1900
#             pointing_table[i]['datetime'] = \
#                 datetime.datetime(century+int(line[6:8]), int(line[3:5]), int(line[0:2]),
#                                   int(line[9:11]), int(line[12:14]), int(line[15:17]))
#             pointing_table[i]['phase_RA'] = []
#             pointing_table[i]['phase_LA'] = []
#             for j in range(0, 8):
#                 pointing_table[i]['phase_RA'].append(int(line[18+j]))
#                 pointing_table[i]['phase_LA'].append(int(line[48+j]))
#             pointing_table[i]['delay_EW_RA'] = int(line[27:30])
#             pointing_table[i]['delay_EW_LA'] = int(line[57:60])
#             pointing_table[i]['delay_NS_RA'] = int(line[31:34])
#             pointing_table[i]['delay_NS_LA'] = int(line[61:64])
#             pointing_table[i]['Field_filter_RA'] = int(line[35])
#             pointing_table[i]['Field_filter_LA'] = int(line[55])
#             pointing_table[i]['Lab_filter_RA'] = []
#             pointing_table[i]['Lab_filter_LA'] = []
#             for j in range(0, 4):
#                 pointing_table[i]['Lab_filter_RA'].append(int(line[37+j*2]))
#                 pointing_table[i]['Lab_filter_LA'].append(int(line[67+j*2]))
#             if line[45:47] == 'OF':
#                 pointing_table[i]['noise_at_RA'] = 255
#             else:
#                 pointing_table[i]['noise_at_RA'] = int(line[45:47])
#             if line[75:77] == 'OF':
#                 pointing_table[i]['noise_at_LA'] = 255
#             else:
#                 pointing_table[i]['noise_at_LA'] = int(line[75:77])
#             pointing_table[i]['Interarray_delay'] = int(line[78:81])
#             i += 1
#
#     cnt_pou = i
#     return pointing_table


# def rt1_to_ddr_calib(file_data, file_pointing, ignore_first_cal):
#
#     edr = build_edr_data(file_data)
#     out = edr
#     # header, date, data, status_code, index_lh, index_rh
#     pointing = load_strpou_data(file_pointing)
#
#     cnt_lh = sum(edr['index_lh'])
#     cnt_rh = sum(edr['index_rh'])
#
#     print 'LH data :'+str(cnt_lh)+' sweeps'
#     print 'RH data :'+str(cnt_rh)+' sweeps'
#
#     print 'Calibration sequence identification...'
#     index_cal = [(x['noise_at_LA'] != 255) | (x['noise_at_RA'] != 255) for x in pointing]
#     cnt_cal = sum([int(i) for i in index_cal])
#     print 'Found '+str(cnt_cal)+' calibration intervals.'
#
#     cal_sequences = dict()
#     cal_sequences['start_time'] = []
#     cal_sequences['stop_time'] = []
#     for x in pointing:
#         if x['noise_at_LA'] == 30:
#             cal_sequences['start_time'].append(dict())
#             cal_sequences['start_time']['']
#         cal_sequences[int(i)]
#
#     if ignore_first_cal:
#         print "IGNORE_FIRST_CAL: not yet implemented"
#
#     return out


