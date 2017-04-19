from spacepy import pycdf
import os
import sys
import struct
import datetime
import numpy as np
import pdb


def edr_cdf_to_rt1(cdf_file, rt1_file='', verbose=False):
    # pdb.set_trace()
    print "Input CDF file: {}".format(os.path.basename(cdf_file))

    # opening CDF
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

    # building datetime object for meridian transit
    meridian_dt = datetime.datetime.strptime(cdf.attrs['NDA_meridian_time'][0], "%Y-%m-%dT%H:%M:%SZ")
    # fixing 2 digits year number
    meridian_yy = meridian_dt.year - 2000
    if meridian_yy < 0:
        meridian_yy += 100

    # building datetime object for observation start
    start_dt = datetime.datetime.strptime(cdf.attrs['PDS_Observation_start_time'][0][0:19], "%Y-%m-%dT%H:%M:%S")
    # fixing 2 digits year number
    start_yy = start_dt.year - 2000
    if start_yy < 0:
        start_yy += 100

    # building datetime object for observation stop
    stop_dt = datetime.datetime.strptime(cdf.attrs['PDS_Observation_stop_time'][0][0:19], "%Y-%m-%dT%H:%M:%S")

    # setting up RF filter header
    n_rf_filter = len(cdf.attrs['NDA_rf_filter_selected'])  # Number of RF filter change in CDF header

    rf0_sele = cdf.attrs['NDA_rf_filter_selected'][0]  # Code of RF filter selected at Observation start
    rf0_dt = datetime.datetime.strptime(cdf.attrs['NDA_rf_filter_time_change'][0][0:19], "%Y-%m-%dT%H:%M:%S")
    rf0_hour = rf0_dt.hour  # Start hour of RF filter 0 (start of observation)
    rf0_minu = rf0_dt.minute  # Start minute of RF filter 0 (start of observation)

    # 1st RF filter change (if exists)
    if n_rf_filter > 1:
        rf1_sele = cdf.attrs['NDA_rf_filter_selected'][1]
        rf1_dt = datetime.datetime.strptime(cdf.attrs['NDA_rf_filter_time_change'][1][0:19], "%Y-%m-%dT%H:%M:%S")
        rf1_hour = rf1_dt.hour
        rf1_minu = rf1_dt.minute
    else:
        rf1_sele = 0
        rf1_hour = 0
        rf1_minu = 0

    # 2nd RF filter change (if exists)
    if n_rf_filter == 3:
        rf2_sele = cdf.attrs['NDA_rf_filter_selected'][2]
        rf2_dt = datetime.datetime.strptime(cdf.attrs['NDA_rf_filter_time_change'][1][0:19], "%Y-%m-%dT%H:%M:%S")
        rf2_hour = rf2_dt.hour
        rf2_minu = rf2_dt.minute
    else:
        rf2_sele = 0
        rf2_hour = 0
        rf2_minu = 0

    # building header (version 6, see srn_nda_routine_jup.py)
    header = 'x{:02d}{:02d}{:03d}{:03d}{:05d}{:02d}{:02d}{:02d}{:1d}{:02d}:{:02d}{:1d}{:02d}:{:02d}{:1d}{:02d}:{:02d}' \
             ' {:02d}/{:02d}/{:02d} {:02d}/{:02d}/{:02d} {:02d}:{:02d}'.format(
        int(np.round(cdf.attrs['VESPA_spectral_range_min'][0] / 1e6)),  # Freq min (MHz)
        int(np.round(cdf.attrs['VESPA_spectral_range_max'][0] / 1e6)),  # Freq max (MHz)
        int(np.round(cdf.attrs['VESPA_spectral_resolution'][0] / 1e3)),  # Spectral Res (kHz)
        int(np.round(cdf.attrs['NDA_reference_level'][0])),  # Ref level (dB)
        int(np.round(cdf.attrs['NDA_sweep_duration'][0] * 1e3)),  # Sweep duration (ms)
        int(np.round(cdf.attrs['NDA_power_resolution'][0])),  # Power resolution (dB/div)
        meridian_dt.hour, meridian_dt.minute,  # Meridian Transit time: Hour, Minute
        rf0_sele, rf0_hour, rf0_minu,  # 1st RF filter selected and time change (hout, minute)
        rf1_sele, rf1_hour, rf1_minu,  # 2nd RF filter selected and time change (hout, minute)
        rf2_sele, rf2_hour, rf2_minu,  # 3rd RF filter selected and time change (hout, minute)
        meridian_dt.day, meridian_dt.month, meridian_yy,  # Meridian Transit time date: Day, Month, Year
        start_dt.day, start_dt.month, start_yy,  # Observation Start date: Day, Month, Year
        stop_dt.hour, stop_dt.minute  # Observation stop time: Hour, Minute
    )
    if verbose:
        print 'Header = {}'.format(header)
        filesize_prev = 0

    # If no RT1 file name provide, building RT1 file name
    if rt1_file == '':
        rt1_file = 'J{:02d}{:02d}{:02d}.RT1'.format(meridian_yy, meridian_dt.month, meridian_dt.day)
    print "Output RT1 file: {}".format(os.path.basename(rt1_file))

    # open RT1 file in write binary mode
    rt1 = open(rt1_file, 'wb')

    # writing header, with trailing null characters
    rt1_header = '{:%<405s}'.format(header).replace('%',chr(0))
    rt1.write(rt1_header)

    data_dt = cdf['Epoch']
    n_data = len(data_dt)

    # looping on sweeps
    for ii in range(n_data):
        if verbose:
             print 'CDF record id = {:06d}'.format(ii)

        # LL sweep
        rec_dt = data_dt[ii]
        rectime = list()
        rectime.append(rec_dt.hour)
        rectime.append(rec_dt.minute)
        rectime.append(rec_dt.second)
        rectime.append(int(rec_dt.microsecond / 1e4))
        rt1.write(bytearray(rectime))
        rt1.write(bytearray(cdf['LL'][ii]))
        rt1.write(chr(cdf['STATUS'][ii][0]))

        # RR sweep
        rec_dt = data_dt[ii] + datetime.timedelta(seconds=float(cdf['RR_SWEEP_TIME_OFFSET'][ii]))
        rectime = list()
        rectime.append(rec_dt.hour)
        rectime.append(rec_dt.minute)
        rectime.append(rec_dt.second)
        rectime.append(int(rec_dt.microsecond / 1e4))
        rt1.write(bytearray(rectime))
        rt1.write(bytearray(cdf['RR'][ii]))
        rt1.write(chr(cdf['STATUS'][ii][1]))

    rt1.close()
    cdf.close()

    return 0


if __name__ == '__main__':

    arg = sys.argv

    if len(arg) >= 2:

        file_cdf = arg[1]

        if not file_cdf.to_lower.endswith('.cdf'):
            print "Wrong input file. File extension not matching CDF file."
            sys.exit(-1)

        if len(arg) == 3:

            file_rt1 = arg[2]

            if not file_rt1.endswidth('.RT1'):
                print "Wrong output file. File extension must be '.RT1'."
                sys.exit(-1)

            edr_cdf_to_rt1(file_cdf, file_rt1)

        else:

            edr_cdf_to_rt1(file_cdf)

        sys.exit(0)

    else:

        print "Wrong syntax."
        print "Usage: edr_cdf_to_rt1 input.cdf [OUTPUT.RT1]"
        sys.exit(-1)
