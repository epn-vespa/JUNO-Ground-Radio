#! /usr/bin/python
from spacepy import pycdf
from astropy.time import Time
import os
import platform
import sys
import psycopg2
import datetime
import hashlib
import json

def init_md():
    """
    Initialize metadata object
    :return md: metadata object
    """

    md = dict()
    md['len'] = 0
    md['tab'] = []
    return md

def md5sum(filename, blocksize=65536):
    """
    compute MD5 hash for file
    :param filename: input filename
    :param blocksize: size of block for split computing
    :return digest: MD5 hash for input file
    """
    md5hash = hashlib.md5()
    with open(filename, "rb") as fhdl:
        for block in iter(lambda: fhdl.read(blocksize), b""):
            md5hash.update(block)
    return md5hash.hexdigest()
    
def creation_date(path_to_file):
    """
    Try to get the date that a file was created, falling back to when it was
    last modified if that isn't possible.
    See http://stackoverflow.com/a/39501288/1709587 for explanation.
    :param path_to_file: path to file
    :return datetime: creation date of file
    """

    if platform.system() == 'Windows':
        t = os.path.getctime(path_to_file)
    else:
        stat = os.stat(path_to_file)
        try:
            t=  stat.st_birthtime
        except AttributeError:
            # We're probably on Linux. No easy way to get creation dates here,
            # so we'll settle for when its content was last modified.
            t = stat.st_mtime
    return datetime.datetime.fromtimestamp(t).isoformat()
    
def modification_date(filename):
    """
    extract modification date
    :param filename: path to file
    :return datetime: modification date
    """

    t = os.path.getmtime(filename)
    return datetime.datetime.fromtimestamp(t).isoformat()
    
def access_url(file_name, parent_name, file_type):
    """
    assemble access URL
    :param file_name: file name (no path)
    :param parent_name: name of parent file
    :param file_type: type of access url to be built
    :return url: access file url
    """
    
    yy = int(parent_name[1:3])
    if yy > 70:
        yy +=1900
    else: 
        yy +=2000
    mm = int(parent_name[3:5])
    if file_type == 'cdf':
        root_url = 'http://realtime.obs-nancay.fr/routine_jupiter'
        return '{}/{:04d}/{:02d}/{}'.format(root_url, yy, mm, file_name)
    elif file_type == 'pdf':
        root_url = 'http://realtime.obs-nancay.fr/ql_routine_jupiter/B-W'
        return '{}/{:04d}/{:02d}/{}'.format(root_url, yy, mm, parent_name+'.pdf')
    elif file_type == 'rt1':
        root_url = 'http://realtime.obs-nancay.fr/routine_jupiter'
        return '{}/{:04d}/{:02d}/{}'.format(root_url, yy, mm, parent_name+'.RT1')
    elif file_type == 'png':
        root_url = 'http://realtime.obs-nancay.fr/ql_routine_jupiter/B-W'
        return '{}/{:04d}/{:02d}/{}'.format(root_url, yy, mm, 'T_'+parent_name+'.png')
    else:
        print 'Illegal file type... Aborting'
        return '' 

def nda_metadata(file_cdf, file_pdf, type):
    """
    extract metadata from CDF
    :param file_cdf: input CDF file path
    :param file_pdf: input PDF file path
    :return metadata: metadata for file
    :return status: 0 if success; -1 if error
    """
    
    cdf = pycdf.CDF(file_cdf)
    metadata = {}
    if type == 'cdf':
        cur_file = file_cdf
        metadata['granule_uid'] = cdf.attrs['Parents'][0].split('.')[0]+'_cdf'
        metadata['granule_gid'] = 'SRN NDA EDR Routine Data'
        metadata['creation_date'] = creation_date(cur_file)
        metadata['access_format'] = cdf.attrs['VESPA_access_format'][0]
        metadata['processing_level'] = 2
    elif type == 'pdf':
        cur_file = file_pdf
        metadata['granule_uid'] = cdf.attrs['Parents'][0].split('.')[0]+'_pdf'
        metadata['granule_gid'] = 'SRN NDA PDF Routine Summary Plots'
        metadata['creation_date'] = cdf.attrs['PDS_Observation_stop_time'][0]
        metadata['access_format'] = 'application/pdf'
        metadata['processing_level'] = 4
    elif type == 'rt1':
        cur_file = os.path.join(os.path.dirname(file_cdf), cdf.attrs['Parents'][0])
        metadata['granule_uid'] = cdf.attrs['Parents'][0].split('.')[0]
        metadata['granule_gid'] = 'SRN NDA RT1 Routine Raw Data'
        metadata['creation_date'] = cdf.attrs['PDS_Observation_stop_time'][0]
        metadata['access_format'] = 'application/binary'
        metadata['processing_level'] = 1
    else:
        print 'Illegal file type... Aborting'
        return metadata,-1
        
    metadata['obs_id'] = cdf.attrs['Parents'][0].split('.')[0]
    metadata['target_name'] = '#'.join(cdf.attrs['PDS_Observation_target'])
    metadata['target_class'] = '#'.join(cdf.attrs['VESPA_target_class'])
    metadata['time_min'] = float(Time(cdf.attrs['PDS_Observation_start_time']).jd[0])
    metadata['time_max'] = float(Time(cdf.attrs['PDS_Observation_stop_time']).jd[0])
    metadata['time_sampling_step'] = float(cdf.attrs['VESPA_time_sampling_step'][0])
    metadata['time_exp'] = float(cdf.attrs['VESPA_time_exp'][0])
    metadata['spectral_range_min'] = float(cdf.attrs['VESPA_spectral_range_min'][0])
    metadata['spectral_range_max'] = float(cdf.attrs['VESPA_spectral_range_max'][0])
    metadata['spectral_sampling_step'] = float(cdf.attrs['VESPA_spectral_sampling_step'][0])
    metadata['spectral_resolution'] = float(cdf.attrs['VESPA_spectral_resolution'][0])
    metadata['measurement_type'] = '#'.join(cdf.attrs['VESPA_measurement_type'])
    metadata['receiver_name'] = cdf.attrs['VESPA_receiver_name'][0].split('>')[0]
    metadata['modification_date'] = modification_date(cur_file)
    metadata['release_date'] = datetime.datetime.now().isoformat()
    metadata['access_md5'] = md5sum(cur_file)
    metadata['thumbnail_url'] = access_url(os.path.basename(file_cdf), data['obs_id'], 'png')
    metadata['target_region'] = '#'.join(cdf.attrs['VESPA_target_region'])
    metadata['feature_name'] = '#'.join(cdf.attrs['VESPA_feature_name'])
    metadata['bib_reference'] = '#'.join(cdf.attrs['VESPA_bib_reference'])
    metadata['access_url'] = access_url(os.path.basename(file_cdf), data['obs_id'], type)
    metadata['file_name'] = os.path.basename(cur_file)
    metadata['access_estsize'] = os.path.getsize(cur_file)/1024
    cdf.close()

    return metadata,0

def ingest(md, cur):
    """
    ingest metadata into nda.metadata
    :param md: metadata of file
    :param cur: current database access handle
    """
    
    for item in md['tab']:
        cur.execute("select count(*) from nda.metadata where granule_uid = {}".format(item['granule_uid']) # test to avoid multiple same entry
        number_of_rows = cur.fetchone()[0]
        if number_of_rows == 1:
            query = "UPDATE nda.metadata SET " 
            first = True
            for k in item.keys(): 
                if first:
                    query += "{}={}".format(k,'%s')
                    param = (item[k],)
                    first = False
                else:
                    query += ", {}={}".format(k,'%s')
                    param = param + (item[k],)
            query += " WHERE granule_uid = %s;"
            param = param + (item['granule_uid'],)
            
        else:
            query = "INSERT INTO nda.metadata "  
            query_key = ""
            query_val = ""
            first = True
            for k in item.keys():
                if first:
                    query_key += "({}".format(k)
                    query_val += "(%s"
                    param = (item[k],)
                    first = False
                else:
                    query_key += ", {}".format(k)
                    query_val += ", %s"
                    param = param + (item[k],)
            query_key += ")"
            query_val += ")"
            query += """{} VALUES {};""".format(query_key,query_val) 

        print cur.mogrify(query, param)
        cur.execute(query, param)

    return

file_cdf = sys.argv[1]
file_pdf = sys.argv[2]

if file_cdf == '--help':
    print "NAME"
    print "    upsert_nda.py  -- update or insert nda record in nda.metadata table"
    print "SYNOPSIS"
    print "    upsert_nda.py cdf_file pdf_file"
    print "DESCRIPTION"
    print "    This script is filling the NDA.METADATA table in vogate.obs-nancay.fr."
    print "    - cdf_file: full path to CDF file record to be updated or inserted."
    print "    - pdf_file: full path to PDF summary plot file corresponding to CDF data file."
    print "DEPENDENCIES"
    print "    The following python libraries must be installed: spacepy, astropy, psycopg2, datetime, hashlib, json"

metadata=init_md()
file_types = ['cdf', 'pdf', 'rt1']

for type in file_types:
    data,status = nda_metadata(file_cdf, file_pdf, type)
    if status == 0:
        metadata['tab'].append(data)
        metadata['len'] += 1

# print metadata

# copy the file credential.json to credential_vogate.json and insert gavoadmin password in file.
with open('credential_vogate.json') as f:
    mydb = json.load(f)

mydatabase = mydb['database']
myhost = mydb['host']
myuser = mydb['user']
mypasswd = mydb['passwd']

try:
    con = psycopg2.connect(database= mydatabase, user= myuser, host= myhost, password= mypasswd, port='5432')
    cur = con.cursor()

except psycopg2.DatabaseError, e:
    print 'Error %s' % e
    sys.exit(1)

ingest(metadata, cur)

if con:
    con.commit()
    con.close()
else:
    print "Problem to access the database, check the configuration"


