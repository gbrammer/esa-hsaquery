"""
Fetch datasets selected with the ESA query with MAST Mashup tools

https://mast.stsci.edu/api/v0/MastApiTutorial.html
"""
import sys
import os
import time
import re
import json
import shutil

try: # Python 3.x
    from urllib.parse import quote as urlencode
    from urllib.request import urlretrieve
except ImportError:  # Python 2.x
    from urllib import pathname2url as urlencode
    from urllib import urlretrieve

try: # Python 3.x
    import http.client as httplib 
except ImportError:  # Python 2.x
    import httplib   

from astropy.table import Table
import numpy as np

import pprint
pp = pprint.PrettyPrinter(indent=4)

### Mashup / MAST
def mastQuery(request):
    """Perform a MAST query.
    
        Parameters
        ----------
        request (dictionary): The MAST request json object
        
        Returns head,content where head is the response HTTP headers, and content is the returned data"""
    
    server='mast.stsci.edu'

    # Grab Python Version 
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain",
               "User-agent":"python-requests/"+version}

    # Encoding the request as a json string
    requestString = json.dumps(request)
    requestString = urlencode(requestString)
    
    # opening the https connection
    conn = httplib.HTTPSConnection(server)

    # Making the query
    conn.request("POST", "/api/v0/invoke", "request="+requestString, headers)

    # Getting the response
    resp = conn.getresponse()
    head = resp.getheaders()
    content = resp.read().decode('utf-8')

    # Close the https connection
    conn.close()

    return head,content

try:
    from .fetch import DEFAULT_PRODUCTS
except:
    from hsaquery.fetch import DEFAULT_PRODUCTS
    
def directDownload(accessLink, filename, path='./'):
    """
    https://mast.stsci.edu/api/v0/pyex.html
    """
    
    server='mast.stsci.edu'
    conn = httplib.HTTPSConnection(server)
    conn.request("GET", "/api/v0/download/file/"+accessLink.lstrip('mast:'))
    resp = conn.getresponse()
    fileContent = resp.read()
    
    with open(os.path.join(path, filename),'wb') as FLE:
        FLE.write(fileContent)
    
    conn.close()
    
def get_from_MAST(table, inst_products=DEFAULT_PRODUCTS, zipFilename='mastDownload', request_only=False, retrieve=True, direct=False, path='./', skip_existing=True):
    """
    testing
    """
    
    import numpy as np
    import os

    ## xx testing
    if False:
        from hsaquery import query, overlaps
        from grizli.pipeline import auto_script
        table = query.run_query(box=None, proposid=[11359], instruments=['WFC3', 'ACS'], extensions=['FLT'], filters=['G102','G141'], extra=[])
        
    # Compute the products to fetch
    products = []
    for i in range(len(table)):
        inst_det = '{0}/{1}'.format(table['instrument'][i], table['detector'][i]) 
        
        if inst_det in inst_products:
            products.append(inst_products[inst_det])
        else:
            products.append(['RAW'])
    
    URLs = []
    productTypes = []
    outPaths = []
    descriptions = []
    
    for i, obs in enumerate(table['observation_id']):
        for p in products[i]:
            URLs.append( 'mast:HST/product/{0}/{0}_{1}.fits'.format(obs.lower(), p.lower()))
            
            productTypes.append('image')
            outPaths.append( '{0}/{1}_{2}.fits'.format(path, obs.lower(), p.lower()))
     
    if direct:
        filenames = [os.path.basename(file) for file in outPaths]
        il = range(len(filenames))
        for i, accessLink, filename in zip(il, URLs, filenames):
            print('({0:>3d}/{1:>3d}): {2}'.format(i+1, il[-1]+1, os.path.join(path, filename)))
            if os.path.exists(os.path.join(path, filename)) & skip_existing:
                continue
                
            directDownload(accessLink, filename, path=path)
        
        return True
        
    #zipFilename = "mastDownload"
    extension = "tar.gz"

    # prepare the bundle
    bundleRequest = {"service":"Mast.Bundle.Request",
                     "params":{"urlList":",".join(URLs),
                               "filename":zipFilename,
                               "pathList":",".join(outPaths),
                               "descriptionList":list(descriptions),
                               "productTypeList":list(productTypes),
                               "extension":extension},
                     "format":"json",
                     "page":1,
                     "pagesize":1000}  

    headers,bundleString = mastQuery(bundleRequest)
    bundleInfo = json.loads(bundleString)

    #pp.pprint(bundleInfo)
    
    if retrieve:
        # Fetch it
        print('Retrieve {0}'.format(bundleInfo['url']))
        urlretrieve(bundleInfo['url'], zipFilename+"."+extension)
    
    print(bundleInfo['url'])
    return(bundleInfo)
        