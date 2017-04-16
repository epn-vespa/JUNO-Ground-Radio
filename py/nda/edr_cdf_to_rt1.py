from spacepy import pycdf
import os, sys
import struct
import datetime
import numpy as np

def edr_cdf_to_rt1(cdf_file, rt1_file=''):

    print,"Input CDF file: {}".format(os.path.basename(cdf_file))

    cdf = pycdf.CDF(cdf_file)

    # header['freq_min'] = raw[1:3]    # MHz
    # header['freq_max'] = raw[3:5]    # MHz
    # header['freq_res'] = raw[5:8]    # kHz
    # header['ref_levl'] = raw[8:11]   # dBm
    # header['swp_time'] = raw[11:16]  # ms
    # header['powr_res'] = raw[16:18]  # dB/div
    # header['merid_hh'] = raw[18:20]  # Meridian time (hours)
    # header['merid_mm'] = raw[20:22]  # Meridian time (minutes)
    # header['rf0_sele'] = raw[22:23]  # RF Filter 0 (selected at start of observations)
    # header['rf0_hour'] = raw[23:25]  # Start for RF filter 0 (hours)
    # header['rf0_minu'] = raw[26:28]  # Start for RF filter 0 (minutes)
    # header['rf1_sele'] = raw[28:29]  # RF Filter 1
    # header['rf1_hour'] = raw[29:31]  # Start for RF filter 1 (hours)
    # header['rf1_minu'] = raw[32:34]  # Start for RF filter 1 (minutes)
    # header['rf2_sele'] = raw[34:35]  # RF Filter 2
    # header['rf2_hour'] = raw[35:37]  # Start for RF filter 2 (hours)
    # header['rf2_minu'] = raw[38:40]  # Start for RF filter 2 (minutes)
    # header['merid_dd'] = raw[41:43]  # Meridian date (day)
    # header['merid_mo'] = raw[44:46]  # Meridian date (month)
    # header['merid_yr'] = raw[47:49]  # Meridian date (year)
    # header['start_dd'] = raw[50:52]  # Observation start date (day)
    # header['start_mo'] = raw[53:55]  # Observation start date (month)
    # header['start_yr'] = raw[56:58]  # Observation start date (year)
    # header['h_stp_hh'] = raw[59:61]  # Observation stop time (hours)
    # header['h_stp_mm'] = raw[62:64]  # Observation stop time (minutes)

    meridian_dt = datetime.datetime.strptime(cdf.attrs['NDA_meridian_time'][0], "%Y-%m-%dT%H:%M:%SZ")

    n_rf_filter = len(cdf.attrs['NDA_rf_filter_selected'])

    rf0_sele = cdf.attrs['NDA_rf_filter_selected'][0]
    rf0_hour = datetime.datetime.strptime(cdf.attrs['NDA_rf_filter_time_change'][0], "%Y-%m-%dT%H:%M:%SZ").hour
    rf0_minu = datetime.datetime.strptime(cdf.attrs['NDA_rf_filter_time_change'][0], "%Y-%m-%dT%H:%M:%SZ").minute

    if n_rf_filter > 1:
        rf1_sele = cdf.attrs['NDA_rf_filter_selected'][1]
        rf1_hour = datetime.datetime.strptime(cdf.attrs['NDA_rf_filter_time_change'][1], "%Y-%m-%dT%H:%M:%SZ").hour
        rf1_minu = datetime.datetime.strptime(cdf.attrs['NDA_rf_filter_time_change'][1], "%Y-%m-%dT%H:%M:%SZ").minute
    else:
        rf1_sele = 0
        rf1_hour = 0
        rf1_minu = 0

    if n_rf_filter == 3:
        rf2_sele = cdf.attrs['NDA_rf_filter_selected'][2]
        rf2_hour = datetime.datetime.strptime(cdf.attrs['NDA_rf_filter_time_change'][2], "%Y-%m-%dT%H:%M:%SZ").hour
        rf2_minu = datetime.datetime.strptime(cdf.attrs['NDA_rf_filter_time_change'][2], "%Y-%m-%dT%H:%M:%SZ").minute
    else:
        rf2_sele = 0
        rf2_hour = 0
        rf2_minu = 0

    header = 'x{:02d}{:02d}{:03d}{:03d}{:05d}{:02d}{:02d}{:02d}{:1d}{:02d}:{:02d}{:1d}{:02d}:{:02d}{:1d}{:02d}:{:02d} {:02d}/{:02d}/{:02d}'.format(
        int(np.round(cdf.attrs['VESPA_spectral_range_min'][0] / 1e6)),
        int(np.round(cdf.attrs['VESPA_spectral_range_max'][0] / 1e6)),
        int(np.round(cdf.attrs['VESPA_spectral_resolution'][0] / 1e3)),
        int(np.round(cdf.attrs['NDA_reference_level'][0])),
        int(np.round(cdf.attrs['NDA_sweep_duration'][0] * 1e3)),
        int(np.round(cdf.attrs['NDA_power_resolution'][0])),
        meridian_dt.hour, meridian_dt.minute,
        rf0_sele, rf0_hour, rf0_minu,
        rf1_sele, rf1_hour, rf1_minu,
        rf2_sele, rf2_hour, rf2_minu,
        meridian_dt.day, meridian_dt.month, meridian_dt.year
    )

    if rt1_file == '':
        rt1_file = ''

    rt1 = open(rt1_file,'wb')
    rt1.write(struct.pack(header))

    return 0

if __name__ == '__main__':

    arg = sys.argv

    if len(arg) >= 2:

        file_cdf = arg[1]

        if not file_cdf.to_lower.endswith('.cdf'):

           print,"Wrong input file. File extension not matching CDF file."
           return -1

        if len(arg) == 3:

            file_rt1 = arg[2]
            if not file_rt1.endswidth('.RT1'):

               print,"Wrong output file. File extension must be '.RT1'."
               return -1

        else:
            edr_cdf_to_rt1(file_cdf)

        edr_cdf_to_rt1(file_cdf, file_rt1)

    else:

        print, "Wrong syntax."
        print, "Usage: edr_cdf_to_rt1 input.cdf [OUTPUT.RT1]"
        return -1

