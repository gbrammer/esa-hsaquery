"""
Fetch data directly from ESA Hubble Science Archive
"""

DEFAULT_PRODUCTS = {'WFC3/IR':['RAW'],
                    'WFPC2/WFPC2':['C0M','C1M'],
                    'ACS/WFC':['FLC'],
                    'WFC3/UVIS':['FLC']}
                    
def make_curl_script(table, level=None, script_name=None, inst_products=DEFAULT_PRODUCTS):
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
                curl_list.append('curl {0}{1}_{2} -o {3}_{4}.fits.gz'.format(BASE_URL, o, product, o.lower(), product.lower()))
            
    else:
        curl_list = ['curl {0}{1}_{2} -o {3}_{4}.fits.gz'.format(BASE_URL, o, level, o.lower(), level.lower()) for o in table['observation_id']]
    
    if script_name is not None:
        fp = open(script_name, 'w')
        fp.write('\n'.join(curl_list))
        fp.close()
    
    return curl_list
        
    