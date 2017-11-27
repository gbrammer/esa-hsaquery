"""
Fetch data directly from ESA Hubble Science Archive
"""

DEFAULT_PRODUCTS = {'WFC3/IR':['RAW'],
                    'WFPC2/WFPC2':['C0M','C1M'],
                    'ACS/WFC':['FLC'],
                    'WFC3/UVIS':['FLC']}
                    
def make_curl_script(table, level=None, script_name=None, inst_products=DEFAULT_PRODUCTS, skip_existing=True, output_path='./'):
    """
    Generate a "curl" script to fetch products from the ESA HSA
    
    Parameters
    ----------
    table : `~astropy.table.Table`
        Table output from `~hsaquery.query` scripts.
        
    level : str
        Specific data type to retrieve (e.g., 'FLT', 'RAW', 'SPT', etc.).  
        If `None`, then retrieve the following:
            'RAW'       for WFC3/IR
            'FLC'       for ACS/WFC and WFC3/UVIS
            'COM'+'C1M' for WFPC2
            
    script_name : str or None
        If a string, then save the curl commands to a file.
    
    Returns
    -------
    curl_list : list
        List of curl commands.
    
    """
    import tempfile
    import glob
    import os
    
    BASE_URL = 'http://archives.esac.esa.int/ehst-sl-server/servlet/data-action?ARTIFACT_ID=' #J6FL25S4Q_RAW'
        
    if level is None:
        # Get RAW for WFC3/IR, FLC for UVIS and ACS
        curl_list = []
        for i in range(len(table)):
            inst_det = '{0}/{1}'.format(table['instrument'][i], table['detector'][i]) 
            
            if inst_det in inst_products:
                products = inst_products[inst_det]
            else:
                products = ['RAW']
        
            o = table['observation_id'][i]
            for product in products:
                if skip_existing:
                    path = '{2}/{0}_{1}.fits*'.format(o.lower(), product.lower(), output_path)
                    if len(glob.glob(path)) > 0:
                        skip = True
                    else:
                        skip = False
                else:
                    skip = False
                    
                if not skip:
                    curl_list.append('curl {0}{1}_{2} -o {5}/{3}_{4}.fits.gz'.format(BASE_URL, o, product, o.lower(), product.lower(), output_path))
            
    else:
        curl_list = ['curl {0}{1}_{2} -o {3}_{4}.fits.gz'.format(BASE_URL, o, level, o.lower(), level.lower()) for o in table['observation_id']]
    
    if script_name is not None:
        fp = open(script_name, 'w')
        fp.write('\n'.join(curl_list))
        fp.close()
    
    return curl_list

def persistence_products(tab):
    import numpy as np
    wfc3 =  tab['instdet'] == 'WFC3/IR'
    progs = np.unique(tab[wfc3]['proposal_id'])
    persist_files = []
    for prog in progs:
        p = wfc3 & (tab['proposal_id'] == prog)
        visits = np.unique(tab['visit'][p])
        print(visits)
        for visit in visits:
            if isinstance(visit, str):
                vst = visit
            else:
                vst = '{0:02d}'.format(visit)
        
            file_i = 'https://archive.stsci.edu/pub/wfc3_persist/{0}/Visit{1}/{0}.Visit{1}.tar.gz'.format(prog, vst)
            persist_files.append(file_i)

    return persist_files
    
    
    