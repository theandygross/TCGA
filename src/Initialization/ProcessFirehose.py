"""
Created on Jul 5, 2013

Process data from firehose_get bulk download. 
1) get firehose_get script and cd into directory to store data
2) Download data
    !bash firehose_get stddata
    !bash firehose_get analyses
    
3) run process_firehose_get(data_path, cancer, date)

4) delete tar files downloaded from firehose_get

@author: agross
"""
import os as os
import tarfile


def checkMD5(f):
    """http://stackoverflow.com/questions/1131220/
       get-md5-hash-of-a-files-without-open-it-in-python"""
    import hashlib
    md5 = hashlib.md5()
    with open(f, 'rb') as f: 
        for chunk in iter(lambda: f.read(8192), b''): 
            md5.update(chunk)
    return md5.hexdigest()


def check_md5_get_files(path):
    files = os.listdir(path)
    files = [f for f in files if f + '.md5' in files]
    for f in files[::2]:
        with open(path + f + '.md5', 'rb') as infile:
            md5 = infile.read()
        assert checkMD5(path + f) == md5.split(' ')[0]
    return files


def unpack(files, path, new_path):
    for f in files:
        t = tarfile.open(path + f)
        try:
            os.makedirs(new_path)
        except:
            pass
        t.extractall(path=new_path)
    for f in os.listdir(new_path):
        if '.' not in f:
            continue
        try:
            os.rename('{}/{}'.format(new_path, f),
                      '{}/{}'.format(new_path, f.split('.')[3]))
        except OSError:
            print path
            print f.split('.')[3]


def clean_file_tree(new_path):
    """
    Cleans up filenames and directory structure to make it a little more
    human readable. Unflattens directory structure (Broad uses __ to 
    flatten directories).   
    """
    for folder in os.listdir(new_path):  # get rid Merge in folder names
        fp = new_path + '/' + folder + '/'
        if folder.startswith('Merge_'):
            folder = folder[6:]
            os.rename(fp, new_path + '/' + folder)
            fp = new_path + '/' + folder        
        if '__' in folder:  # un-flatten directories
            os.renames(fp, new_path + '/' + '/'.join(folder.split('__')))
            
    for ll in os.walk(new_path):  # get rid of directory names in filenames
        for f in ll[2]:
            if ('__' in f) and (f.startswith('gdac') is False):
                os.rename(ll[0] + '/' + f,
                          ll[0] + '/' + '.'.join(f.split('.')[-2:]))


def process_files(path, new_path):
    """
    Read in files from firehose_get output and convert to 
    more convenient directory structure for data-analysis.
    """
    files = check_md5_get_files(path)
    ffpe = [f for f in files if 'FFPE' in f]
    if len(ffpe) > 0:
        process_files_sub(ffpe, path, new_path + '_FFPE')

    frozen = [f for f in files if 'FFPE' not in f]
    if len(frozen) > 0:
        process_files_sub(frozen, path, new_path)


def process_files_sub(files, path, new_path):
    mage_tabs = [f for f in files if 'mage-tab' in f]
    unpack(mage_tabs, path, '{}/mage-tabs'.format(new_path))

    aux_files = [f for f in files if '.aux.' in f]
    unpack(aux_files, path, '{}/aux-files'.format(new_path))

    reg_files = [f for f in files if ('.aux.' not in f)
                 and ('mage-tab' not in f)]

    unpack(reg_files, path, new_path)


def process_firehose_get(data_path, cancer, date):
    """
    Wrapper to format data from firehose_get download into directory 
    structure for the rest of the analysis. 

    Input:
        data_path: base path where data is stored
        date: date of versioned firehose run in YYYY_MM_DD format.
        cancer: name of cancer to process
    """
    cancer_dir = '{}/analyses__{}/{}/'.format(data_path, date, cancer)
    path = cancer_dir + date.replace('_', '') + '/'
    new_path = '{}/Firehose__{}/analyses/{}'.format(data_path, date, cancer)
    process_files(path, new_path)
    
    cancer_dir = '{}/stddata__{}/{}/'.format(data_path, date, cancer)
    path = cancer_dir + date.replace('_', '') + '/'
    new_path = '{}/Firehose__{}/stddata/{}'.format(data_path, date, cancer)
    process_files(path, new_path)
    clean_file_tree(new_path)


def process_all_cancers(firehose_path, date):
    """
    Process data retrieved using firehose_get.
    
    firehose_path: Path to top level directory where firehose_get was run.
    date: date of versioned firehose run in YYYY_MM_DD format.
    """
    for cancer in os.listdir('{}/analyses__{}'.format(firehose_path, date)):
        if '.' in cancer:  # random files stuck in the directory
            break
        process_firehose_get(firehose_path, cancer, date)
