#!/usr/bin/env python

"""SRN NDA NewRoutine Jup EDR dataset building scripts"""

__author__ = "Baptiste Cecconi"
__copyright__ = "Copyright 2017, LESIA-PADC-USN, Observatoire de Paris"
__credits__ = ["Baptiste Cecconi", "Andree Coffre", "Laurent Lamy"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Baptiste Cecconi"
__email__ = "baptiste.cecconi@obspm.fr"
__status__ = "Development"

import numpy as np
import os
import os.path
import struct
import datetime
from spacepy import pycdf
import json
import pdb


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


def load_dat_data(input_file, nrec=-1):
    """
    Loads RT1 data
    :param input_file: input RT1 file
    :param nrec: number of record to load, default (=-1) means load all records.
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

    if nrec == -1:
        packet_size = 10000
    else:
        packet_size = nrec
        nspectra = nrec

    ifrq_min = 0
    ifrq_max = 0
    for i in range(2048):
        if header['frq'][i] < 10:
            ifrq_min = i + 1
        if header['frq'][i] <= 40:
            ifrq_max = i

    header['freq_min'] = header['frq'][ifrq_min]  # MHz
    header['freq_max'] = header['frq'][ifrq_max]  # MHz
    header['freq_res'] = 48.828125  # kHz
    header['frequncy'] = np.array(header['frq'])[ifrq_min:ifrq_max + 1]
    header['freq_stp'] = 48.828125  # kHz

    data = list()
    jdatetime = list()

    for ii in range(0, nspectra, packet_size):
        read_size = packet_size
        if (ii+packet_size)*record_size+header_size > file_size:
            read_size = nspectra-ii
        print '   {:05d} to {:05d}: reading {} records'.format(ii, ii+read_size, read_size)
        packet_raw = ff.read(record_size*read_size)
        for jj in range(0, read_size):
            rec_val = struct.unpack(record_fmt, packet_raw[jj*record_size:(jj+1)*record_size])
            record = dict()
            record['header'] = dict()
            record['header']['magic'] = rec_val[0]
            magic_code_record = '0x{:08X}'.format(record['header']['magic'])
            if magic_code_record != '0x7F800000':
                print('ERROR at {}'.format(jj))
            record['header']['cnt'] = rec_val[1]
            record['header']['jdatetime'] = rec_val[2:6]
            record['data'] = list()
            for ichan in range(0, nbchan):
                record['data'].append(dict())
                record['data'][ichan]['magic'] = rec_val[8+ichan*2050]
                magic_code_spect = '0x{:08X}'.format(record['data'][ichan]['magic'])
                if magic_code_spect != '0xFF800001':
                    print('ERROR at {}.{}'.format(jj, ichan))
                record['data'][ichan]['spec_typ'] = rec_val[9+ichan*2050]
                record['data'][ichan]['spectrum'] = rec_val[10+ichan*2050+ifrq_min:10+ichan*2050+ifrq_max+1]
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

    # pdb.set_trace()

    var = dict()
    var['data'] = data
    var['jdatetime'] = jdatetime

    header['obsty_id'] = 'srn'
    header['instr_id'] = 'nda'
    header['recvr_id'] = 'newroutine_jup'

    header['time_stp'] = np.median([(var['jdatetime'][i + 1].to_datetime() -
                                     var['jdatetime'][i].to_datetime()).total_seconds()
                                    for i in range(len(var['jdatetime']) - 1)])  # sec
    header['time_exp'] = header['time_stp'] * 1000  # ms

    return header, var


def build_edr_data(input_file):
    """
    Builds EDR data file from RT2 file
    :param input_file: input RT1 file
    :return:
    """

    print 'Loading data...'
    header, var = load_dat_data(input_file)

    nspectra = len(var['jdatetime'])

    out = dict()
    out['header'] = header
    out['time'] = list()
    out['ll_flux'] = list()
    out['rr_flux'] = list()
    out['lr_flux'] = list()
    out['rl_flux'] = list()

    for ii in range(nspectra):
        out['time'].append(var['jdatetime'][ii].to_datetime())
        out['ll_flux'].append(var['data'][ii]['data'][0]['spectrum'])
        out['rr_flux'].append(var['data'][ii]['data'][3]['spectrum'])
        out['lr_flux'].append(var['data'][ii]['data'][1]['spectrum'])
        out['rl_flux'].append(var['data'][ii]['data'][2]['spectrum'])

    return out


def build_edr_cdf(input_file, config_file='config.json', build_cdf_master=False):
    """
    Builds the EDR level CDF file
    :param input_file: input RT1 file
    :param config_file: configuration file (default = config.json in same directory)
    :param build_cdf_master: option to (re)build CDF master (default = False)
    :return:
    """

    cdf_version, sft_version, dat_version = get_versions(config_file)
    master_path, cdfbin_path, cdfout_path = get_paths(config_file, cdf_version)

    edr = build_edr_data(input_file)

    master_cdf = master_cdf_name(edr['header']['obsty_id'], edr['header']['instr_id'], edr['header']['recvr_id'],
                                 'edr', cdf_version)
    skelet_cdf = master_skt_name(edr['header']['obsty_id'], edr['header']['instr_id'], edr['header']['recvr_id'],
                                 'edr', cdf_version)

    master_cdf_path = os.path.join(master_path, master_cdf)
    skelet_cdf_path = os.path.join(master_path, skelet_cdf)

    if build_cdf_master:
        if os.path.exists(master_cdf_path):
            os.remove(master_cdf_path)
        os.system('{} -cdf {} {}'.format(os.path.join(cdfbin_path, 'skeletoncdf'),
                                         master_cdf_path, skelet_cdf_path))

    print "Raw data file name:"
    print input_file
    print "Master CDF file name:"
    print master_cdf

    ndata = len(edr['time'])
    frequency = edr['header']['frequncy']

    cdfout_file = "{}_{}_{}_edr_{:%Y%m%d%H%M}_{:%Y%m%d%H%M}_V{}.cdf"\
        .format(edr['header']['obsty_id'], edr['header']['instr_id'], edr['header']['recvr_id'],
                edr['time'][0], edr['time'][ndata-1], cdf_version)
    cdfout_file_path = os.path.join(cdfout_path, cdfout_file)

    if os.path.exists(cdfout_file_path):
        os.remove(cdfout_file_path)

    print "Output CDF file name:"
    print cdfout_file

    cdfout = pycdf.CDF(cdfout_file_path, master_cdf_path)

    # SETTING PDS GLOBAL ATTRIBUTES
    cdfout.attrs['PDS_Observation_start_time'] = edr['time'][0].isoformat()+'Z'
    cdfout.attrs['PDS_Observation_stop_time'] = edr['time'][ndata-1].isoformat()+'Z'

    # SETTING VESPA GLOBAL ATTRIBUTES
    cdfout.attrs['VESPA_time_sampling_step'] = edr['header']['time_stp']        # in seconds
    cdfout.attrs['VESPA_time_exp'] = edr['header']['time_exp']/1.e3			    # in seconds

    cdfout.attrs['VESPA_spectral_range_min'] = edr['header']['freq_min']*1.e6   # In Hz
    cdfout.attrs['VESPA_spectral_range_max'] = edr['header']['freq_max']*1.e6   # In Hz
    cdfout.attrs['VESPA_spectral_sampling_step'] = edr['header']['freq_stp']*1.e3  # In Hz
    cdfout.attrs['VESPA_spectral_resolution'] = edr['header']['freq_res']*1.e3     # In Hz

    # SETTING SRN-NDA-NewRoutine GLOBAL ATTRIBUTES
    cdfout.attrs['NDA_newroutine_sel_prod'] = [edr['header']['sel_prod0'], edr['header']['sel_prod1']]
    cdfout.attrs['NDA_newroutine_acc_fact'] = edr['header']['acc']

    # SETTING OTHER GLOBAL ATTRIBUTES
    cdfout.attrs['Logical_file_id'] = cdfout_file
    cdfout.attrs['Data_version'] = dat_version
    cdfout.attrs['Skeleton_version'] = cdf_version
    cdfout.attrs['Software_version'] = sft_version
    cdfout.attrs['Software_language'] = 'python'
    cdfout.attrs['Parents'] = os.path.basename(input_file)

    # SETTING VARIABLES
    cdfout['Epoch'] = edr['time']
    cdfout['RR'] = edr['rr_flux']
    cdfout['LL'] = edr['ll_flux']
    cdfout['RL'] = edr['rl_flux']
    cdfout['LR'] = edr['lr_flux']
    cdfout['Frequency'] = frequency

    date_start = edr['time'][0]
    date_stop = edr['time'][ndata-1]
    date_start_round = edr['time'][0].replace(minute=0, second=0, microsecond=0)
    date_stop_round = edr['time'][ndata-1].replace(minute=0, second=0, microsecond=0)+datetime.timedelta(hours=1)

    # SETTING VARIABLES ATTRIBUTES
    cdfout['Epoch'].attrs['VALIDMIN'] = date_start
    cdfout['Epoch'].attrs['VALIDMAX'] = date_stop
    cdfout['Epoch'].attrs['SCALEMIN'] = date_start_round
    cdfout['Epoch'].attrs['SCALEMAX'] = date_stop_round

    cdfout['Frequency'].attrs['SCALEMIN'] = 10.
    cdfout['Frequency'].attrs['SCALEMAX'] = 40.
    cdfout['Frequency'].attrs['UNITS'] = 'MHz'

    cdfout.close()

    print "CDF file written:"
    print cdfout_file_path


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


