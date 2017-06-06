import numpy as np
import os
import os.path
import datetime
import struct as struct
from astropy.time import Time
from spacepy import pycdf


def load_jds_header_type1(hdr_raw):
    hdr_fmt = '<32s32s32s32s8h16b96s256s32i5if24i16s16s32s50i'
    hdr_values = struct.unpack(hdr_fmt,hdr_raw)

    header = dict()
    header['name'] = hdr_values[0].strip('\x00')
    header['time'] = hdr_values[1].strip('\x00')
    header['gmtt'] = hdr_values[2].strip('\x00')
    header['sysn'] = hdr_values[3].strip('\x00')
    header['utc'] = hdr_values[4:12]
    header['syst'] = hdr_values[12:28]
    header['place'] = hdr_values[28].strip('\x00')
    header['desc'] = hdr_values[29].strip('\x00')
    header['PP'] = dict()
    header['PP']['mode'] = hdr_values[30]
    header['PP']['size'] = hdr_values[31]
    header['PP']['Prc_mode'] = hdr_values[32]
    header['PP']['tst_gen'] = hdr_values[33]
    header['PP']['clc'] = hdr_values[34]
    header['PP']['fft_size'] = hdr_values[35]
    header['DSPP'] = dict()
    header['DSPP']['FFT_Size'] = hdr_values[46+16]
    header['DSPP']['MinDSPSize'] = hdr_values[47+16]
    header['DSPP']['MinDMASize'] = hdr_values[48+16]
    header['DSPP']['DMASizeCnt'] = hdr_values[49+16]
    header['DSPP']['DMASize'] = hdr_values[50+16]
    header['DSPP']['Clock_frq'] = hdr_values[51+16]
    header['DSPP']['GPS_Synch'] = hdr_values[52+16]
    header['DSPP']['SSht'] = hdr_values[53+16]
    header['DSPP']['Mode'] = hdr_values[54+16]
    header['DSPP']['WFRch'] = hdr_values[55+16]
    header['DSPP']['Smd'] = hdr_values[56+16]
    header['DSPP']['Offt'] = hdr_values[57+16]
    header['DSPP']['Lb'] = hdr_values[58+16]
    header['DSPP']['Hb'] = hdr_values[59+16]-1
    header['DSPP']['Wb'] = hdr_values[60+16]
    header['DSPP']['NAvr'] = hdr_values[61+16]
    header['DSPP']['CAvr'] = hdr_values[62+16]
    header['DSPP']['Weight'] = hdr_values[63+16]
    header['DSPP']['DCRem'] = hdr_values[64+16]
    header['DSPP']['ExtSyn'] = hdr_values[65+16]
    if header['DSPP']['ExtSyn']:
        header['DSPP']['Clock_frq'] = 10000.*round(header['DSPP']['Clock_frq']/10000.,0)
    header['DSPP']['Ch1'] = hdr_values[66+16]
    header['DSPP']['Ch2'] = hdr_values[67+16]
    header['DSPP']['ExtWin'] = hdr_values[68+16]
    header['DSPP']['Clip'] = hdr_values[69+16]
    header['DSPP']['HPF0'] = hdr_values[70+16]
    header['DSPP']['HPF1'] = hdr_values[71+16]
    header['DSPP']['LPF0'] = hdr_values[72+16]
    header['DSPP']['LPF1'] = hdr_values[73+16]
    header['DSPP']['ATT0'] = hdr_values[74+16]
    header['DSPP']['ATT1'] = hdr_values[75+16]
    header['DSPP']['Soft'] = hdr_values[76+16].strip('\x00')
    header['DSPP']['SVer'] = hdr_values[77+16].strip('\x00')
    header['DSPP']['DSPv'] = hdr_values[78+16].strip('\x00')
    
    header['obsty_id'] = 'utr2'
    header['instr_id'] = 'dspz'
    
    return header


def load_jds_header_type2(hdr_raw):
    hdr_fmt = '<32s32s32s32s8h16b96s256s16i5if24i16s16s32s66i'
    hdr_values = struct.unpack(hdr_fmt,hdr_raw)
    
    header = dict()
    header['name'] = hdr_values[0].strip('\x00')
    header['time'] = hdr_values[1].strip('\x00')
    header['gmtt'] = hdr_values[2].strip('\x00')
    header['sysn'] = hdr_values[3].strip('\x00')
    header['utc'] = hdr_values[4:12]
    header['syst'] = hdr_values[12:28]
    header['place'] = hdr_values[28].strip('\x00')
    header['desc'] = hdr_values[29].strip('\x00')
    header['PP'] = dict()
    header['PP']['mode'] = hdr_values[30]
    header['PP']['size'] = hdr_values[31]
    header['PP']['Prc_mode'] = hdr_values[32]
    header['PP']['tst_gen'] = hdr_values[33]
    header['PP']['clc'] = hdr_values[34]
    header['PP']['fft_size'] = hdr_values[35]
    header['DSPP'] = dict()
    header['DSPP']['FFT_Size'] = hdr_values[46]
    header['DSPP']['MinDSPSize'] = hdr_values[47]
    header['DSPP']['MinDMASize'] = hdr_values[48]
    header['DSPP']['DMASizeCnt'] = hdr_values[49]
    header['DSPP']['DMASize'] = hdr_values[50]
    header['DSPP']['Clock_frq'] = hdr_values[51]
    header['DSPP']['GPS_Synch'] = hdr_values[52]
    header['DSPP']['SSht'] = hdr_values[53]
    header['DSPP']['Mode'] = hdr_values[54]
    header['DSPP']['WFRch'] = hdr_values[55]
    header['DSPP']['Smd'] = hdr_values[56]
    header['DSPP']['Offt'] = hdr_values[57]
    header['DSPP']['Lb'] = hdr_values[58]
    header['DSPP']['Hb'] = hdr_values[59]-1
    header['DSPP']['Wb'] = hdr_values[60]
    header['DSPP']['NAvr'] = hdr_values[61]
    header['DSPP']['CAvr'] = hdr_values[62]
    header['DSPP']['Weight'] = hdr_values[63]
    header['DSPP']['DCRem'] = hdr_values[64]
    header['DSPP']['ExtSyn'] = hdr_values[65]
    if header['DSPP']['ExtSyn']:
        header['DSPP']['Clock_frq'] = 10000.*round(header['DSPP']['Clock_frq']/10000.,0)
    header['DSPP']['Ch1'] = hdr_values[66]
    header['DSPP']['Ch2'] = hdr_values[67]
    header['DSPP']['ExtWin'] = hdr_values[68]
    header['DSPP']['Clip'] = hdr_values[69]
    header['DSPP']['HPF0'] = hdr_values[70]
    header['DSPP']['HPF1'] = hdr_values[71]
    header['DSPP']['LPF0'] = hdr_values[72]
    header['DSPP']['LPF1'] = hdr_values[73]
    header['DSPP']['ATT0'] = hdr_values[74]
    header['DSPP']['ATT1'] = hdr_values[75]
    header['DSPP']['Soft'] = hdr_values[76].strip('\x00')
    header['DSPP']['SVer'] = hdr_values[77].strip('\x00')
    header['DSPP']['DSPv'] = hdr_values[78].strip('\x00')
    
    header['obsty_id'] = 'utr2'
    header['instr_id'] = 'dspz'

    return header


def load_jds_data(file_jds):
#    file_jds = 'C250214_154735.jds'

#     This script was only tested on the file mentioned above ! :-) 
#    waveform and correlation mode are not yet implemented
  
    print 'Loading JDS file [Kharkov/UTR-2/DPS-Z]: '+file_jds
    header_size = 1024
    bps = 4
    bpw = 2
    
#    Getting file size     
    file_size = os.path.getsize(file_jds)

#    Opening file size     
    lun = open(file_jds,'rb')
    raw_header = lun.read(header_size)

#    Loading Header
    if file_jds[0] == 'A' or file_jds[0] == 'B' or file_jds[0] == 'C' or file_jds[0] == 'D' or file_jds[0] == 'E':
        header = load_jds_header_type2(raw_header)
    else: 
        header = load_jds_header_type1(raw_header)

#    Extracting useful variables from Header
    sps = header['DSPP']['Clock_frq']                # Sample per second
    nf = [8192,4096,4096,header['DSPP']['Wb']]       # Number of frequency options
    nfmin0 = [0,0,4096,header['DSPP']['Lb']]         # First frequency options
    nf = nf[header['DSPP']['Offt']]                  # Selected Nf
    nfmin0 = Nfmin0[header['DSPP']['Offt']]          # Selected NFmin0
    nfmax = Nf-1                                     # Last frequency (removing
    nint = header['DSPP']['NAvr']                    # Number of averaged spectrum
    nc = 2-header['DSPP']['Ch1']-header['DSPP']['Ch2']
    
    if header['DSPP']['Mode'] == 0:
        buf_size = 8192*Nc*bpw
    if header['DSPP']['Mode'] == 1:
        buf_size = Nf*Nc*bps      # bytes between 2 consecutive elementary spectra
    if header['DSPP']['Mode'] == 2:
        buf_size = Nf*4*bps
    dt=float(8192)/sps*Nint         # sec / spectrum
    df=sps/float(16384L)       # Hz
    freq=(np.arange(Nf)+Nfmin0)*df/1.e6                   # Frequency ramp (MHz)
    nsweep=float(file_size-header_size)/buf_size
    time=(np.arange(nsweep))*dt                           # Time ramp (sec)
    
#    Loading the rest of the file
    raw_data = lun.read(file_size-header_size)
    lun.close()

#    Conversion of floating point format 
    tmp_data = np.frombuffer(raw_data,dtype='<i').reshape(nsweep,Nf,2)
    
    tmp_data_mant = tmp_data & 0xFFFFFFC0
    tmp_data_sign = 1-(tmp_data & 0x20)/0x10
    tmp_data_expn = (tmp_data & 0x1F)*tmp_data_sign
    tmp_data_snrm = 4*2*1024.0/4294967296.0/Nint
    
    data_float = tmp_data_mant / pow(2.,tmp_data_expn) / tmp_data_snrm
    
#    Creating output variable
    data = {}
    data['spectrum_A']  = data_float[:,0:-1,0]
    data['spectrum_B']  = data_float[:,0:-1,1]
    data['overflow_A']  = (tmp_data[:,-1,0] & 0x80000000)/0x80000000
    data['overflow_B']  = (tmp_data[:,-1,1] & 0x80000000)/0x80000000
    data['microsecond'] = tmp_data[:,-1,0] & 0x7FFFFFFF
    data['fft_count']   = tmp_data[:,-1,1] & 0x7FFFFFFF
            
    return data,freq,time,header
    
def build_edr_data(file):
    
#    Loading EDR data from JDS file
    print 'Loading data...'
    data,freq,time,header = load_jds_data(file)
    ntime = len(time)
    nfreq = len(freq)
    
#    Computing date-times 
    base_date = datetime.datetime(header['utc'][0],header['utc'][1],header['utc'][3],header['utc'][4],header['utc'][5],header['utc'][6],header['utc'][7]*1000)
    data_time = np.array([base_date+datetime.timedelta(seconds=time[i]) for i in range(ntime)])

#    Computing output variable
    out = {}
    out['header'] = header
    out['freq'] = freq
    out['time'] = data_time
    out['data'] = data
    return out

def master_cdf_name(obsty_id, instr_id, level, cdf_version):
    return '{}_{}_{}_000000000000_000000000000_V{}.cdf'.format(obsty_id, instr_id, level, cdf_version)

def master_skt_name(obsty_id, instr_id, level, cdf_version):
    return '{}_{}_{}_000000000000_000000000000_V{}.skt'.format(obsty_id, instr_id, level, cdf_version)

def edr_to_cdf(edr,build_cdf_master):

#    setting up variables
    master_path = 'master/'
    cdfbin_path = '/Applications/cdf/cdf36_0-dist/bin/'
    cdfout_path = 'data/'
    cdf_version = '00'
    dat_version = '00'
    sft_version = '00'

#    Setting SKT and CDF names 
    master_cdf = master_cdf_name(edr['header']['obsty_id'], edr['header']['instr_id'], 'edr', cdf_version)
    skelet_cdf = master_skt_name(edr['header']['obsty_id'], edr['header']['instr_id'], 'edr', cdf_version)
    
#    Creating CDF Master (removing it if already there)
    if build_cdf_master:
        if os.path.exists(master_path+master_cdf):
            os.remove(master_path+master_cdf)
        os.system(cdfbin_path+'skeletoncdf -cdf '+master_path+master_cdf+' '+master_path+skelet_cdf)

    print "Master CDF file name:"
    print master_cdf

#    Creating Time and Frequency Axes 
    ndata = len(edr['time'])
    jul_date = Time(edr['time'],format="datetime",scale="utc").jd.tolist()
    frequency = edr['freq'][0:-1]
    
#    Setting CDF output name 
    cdfout_file = "{}_{}_edr_{:%Y%m%d%H%M}_{:%Y%m%d%H%M}_V{}.cdf".format(edr['header']['obsty_id'], edr['header']['instr_id'], edr['time'][0], edr['time'][ndata-1],cdf_version)
    if os.path.exists(cdfout_path+cdfout_file):
        os.remove(cdfout_path+cdfout_file)

#    Opening CDF object 
    cdfout = pycdf.CDF(cdfout_path+cdfout_file, master_path+master_cdf)

    # SETTING PDS GLOBAL ATTRIBUTES
    cdfout.attrs['PDS_Observation_start_time'] = edr['time'][0].isoformat()+'Z'
    cdfout.attrs['PDS_Observation_stop_time'] = edr['time'][ndata-1].isoformat()+'Z'
    
    # SETTING VESPA GLOBAL ATTRIBUTES
    cdfout.attrs['VESPA_time_min'] = jul_date[0]
    cdfout.attrs['VESPA_time_max'] = jul_date[ndata-1]
    cdfout.attrs['VESPA_time_sampling_step_min'] = np.amin([jul_date[i+1]-jul_date[i] for i in range(0,ndata-2)])*86400.
    cdfout.attrs['VESPA_time_sampling_step_max'] = np.amax([jul_date[i+1]-jul_date[i] for i in range(0,ndata-2)])*86400.
    
    cdfout.attrs['VESPA_spectral_range_min']  = np.amin(frequency)
    cdfout.attrs['VESPA_spectral_range_max']  = np.amax(frequency)
#    cdfout.attrs['VESPA_spectral_sampling_step_min'] = np.amin(frequency[i+1]-frequency[i] for i in range(len(frequency)-1))
#    cdfout.attrs['VESPA_spectral_sampling_step_max'] = np.amax(frequency[i+1]-frequency[i] for i in range(len(frequency)-1))
#    cdfout.attrs['VESPA_spectral_resolution_min'] = 30.e3
#    cdfout.attrs['VESPA_spectral_resolution_max'] = 30.e3
    
    # SETTING OTHER GLOBAL ATTRIBUTES
    cdfout.attrs['Logical_file_id'] = cdfout_file
    cdfout.attrs['Data_version'] = dat_version
    cdfout.attrs['Skeleton_version'] = cdf_version
    cdfout.attrs['Software_version'] = sft_version
    cdfout.attrs['Software_language'] = 'python'
    
    # SETTING UTR-2 GLOBAL ATTRIBUTES
    cdfout.attrs['UTR2_place']      = edr['header']['place']
    cdfout.attrs['UTR2_desc']       = edr['header']['desc']
    cdfout.attrs['UTR2_FFT_Size']   = edr['header']['DSPP']['FFT_Size']   
    cdfout.attrs['UTR2_MinDSPSize'] = edr['header']['DSPP']['MinDSPSize'] 
    cdfout.attrs['UTR2_MinDMASize'] = edr['header']['DSPP']['MinDMASize'] 
    cdfout.attrs['UTR2_DMASizeCnt'] = edr['header']['DSPP']['DMASizeCnt'] 
    cdfout.attrs['UTR2_DMASize']    = edr['header']['DSPP']['DMASize']    
    cdfout.attrs['UTR2_Clock_frq']  = edr['header']['DSPP']['Clock_frq']  
    cdfout.attrs['UTR2_GPS_Synch']  = edr['header']['DSPP']['GPS_Synch']  
    cdfout.attrs['UTR2_SSht']       = edr['header']['DSPP']['SSht']       
    cdfout.attrs['UTR2_Mode']       = edr['header']['DSPP']['Mode']       
    cdfout.attrs['UTR2_WFRch']      = edr['header']['DSPP']['WFRch']      
    cdfout.attrs['UTR2_Smd']        = edr['header']['DSPP']['Smd']        
    cdfout.attrs['UTR2_Offt']       = edr['header']['DSPP']['Offt']       
    cdfout.attrs['UTR2_Lb']         = edr['header']['DSPP']['Lb']         
    cdfout.attrs['UTR2_Hb']         = edr['header']['DSPP']['Hb']         
    cdfout.attrs['UTR2_Wb']         = edr['header']['DSPP']['Wb']         
    cdfout.attrs['UTR2_NAvr']       = edr['header']['DSPP']['NAvr']       
    cdfout.attrs['UTR2_CAvr']       = edr['header']['DSPP']['CAvr']       
    cdfout.attrs['UTR2_Weight']     = edr['header']['DSPP']['Weight']     
    cdfout.attrs['UTR2_DCRem']      = edr['header']['DSPP']['DCRem']      
    cdfout.attrs['UTR2_ExtSyn']     = edr['header']['DSPP']['ExtSyn']     
    cdfout.attrs['UTR2_Ch1']        = edr['header']['DSPP']['Ch1']        
    cdfout.attrs['UTR2_Ch2']        = edr['header']['DSPP']['Ch2']        
    cdfout.attrs['UTR2_ExtWin']     = edr['header']['DSPP']['ExtWin']     
    cdfout.attrs['UTR2_Clip']       = edr['header']['DSPP']['Clip']       
    cdfout.attrs['UTR2_HPF0']       = edr['header']['DSPP']['HPF0']       
    cdfout.attrs['UTR2_HPF1']       = edr['header']['DSPP']['HPF1']       
    cdfout.attrs['UTR2_LPF0']       = edr['header']['DSPP']['LPF0']       
    cdfout.attrs['UTR2_LPF1']       = edr['header']['DSPP']['LPF1']       
    cdfout.attrs['UTR2_ATT0']       = edr['header']['DSPP']['ATT0']       
    cdfout.attrs['UTR2_ATT1']       = edr['header']['DSPP']['ATT1']       
    cdfout.attrs['UTR2_Soft']       = edr['header']['DSPP']['Soft']       
    cdfout.attrs['UTR2_SVer']       = edr['header']['DSPP']['SVer']       
    cdfout.attrs['UTR2_DSPv']       = edr['header']['DSPP']['DSPv']       


    # SETTING VARIABLES
    cdfout['Epoch'] = edr['time']
    cdfout['ISO_DATE'] = [d.isoformat()+'Z' for d in edr['time']]
    cdfout['JD_TIME'] = jul_date
    cdfout['FLUX_A'] = edr['data']['spectrum_A']
    cdfout['FLUX_B'] = edr['data']['spectrum_B']
    cdfout['Frequency'] = frequency*1e6
    
    date_start = edr['time'][0]
    date_stop = edr['time'][ndata-1]
    date_start_round = edr['time'][0].replace(minute=0, second=0, microsecond=0)
    date_stop_round = edr['time'][ndata-1].replace(minute=0, second=0, microsecond=0)+datetime.timedelta(hours=1)
    
    # SETTING VARIABLES ATTRIBUTES
    cdfout['Epoch'].attrs['VALIDMIN'] = date_start
    cdfout['Epoch'].attrs['VALIDMAX'] = date_stop
    cdfout['Epoch'].attrs['SCALEMIN'] = date_start_round
    cdfout['Epoch'].attrs['SCALEMAX'] = date_stop_round
    
    cdfout['ISO_DATE'].attrs['VALIDMIN'] = date_start.isoformat()+'Z'
    cdfout['ISO_DATE'].attrs['VALIDMAX'] = date_stop.isoformat()+'Z'
    cdfout['ISO_DATE'].attrs['SCALEMIN'] = date_start_round.isoformat()+'Z'
    cdfout['ISO_DATE'].attrs['SCALEMAX'] = date_stop_round.isoformat()+'Z'
    
    cdfout['JD_TIME'].attrs['VALIDMIN'] = Time(date_start,format="datetime",scale="utc").jd
    cdfout['JD_TIME'].attrs['VALIDMAX'] = Time(date_stop,format="datetime",scale="utc").jd
    cdfout['JD_TIME'].attrs['SCALEMIN'] = Time(date_start_round,format="datetime",scale="utc").jd
    cdfout['JD_TIME'].attrs['SCALEMAX'] = Time(date_stop_round,format="datetime",scale="utc").jd
    
    cdfout.close()

    
    
