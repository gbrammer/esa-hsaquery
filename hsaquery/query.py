"""
Demo of querying the ESA HST archive
"""

import numpy as np

MASTER_COLORS = {'G102':'#1f77b4',
'F125W':'#ff7f0e',
'F160W':'#2ca02c',
'G141':'#d62728',
'F140W':'#9467bd',
'F105W':'#8c564b',
'F775W':'#8c564b'}

DEFAULT_FIELDS = """OBSERVATION
TARGET
POSITION.RA
POSITION.DEC
POSITION.ECL_LAT
POSITION.ECL_LON
POSITION.GAL_LAT
POSITION.GAL_LON
POSITION.STC_S
POSITION.FOV_SIZE
POSITION.SPATIAL_RESOLUTION
INSTRUMENT.INSTRUMENT_NAME
DETECTOR.DETECTOR_NAME
OPTICAL_ELEMENT.OPTICAL_ELEMENT_NAME
PROPOSAL.PROPOSAL_ID
PROPOSAL.SCIENCE_CATEGORY
PROPOSAL.PI_NAME
ARTIFACT.FILE_NAME
ARTIFACT.FILE_EXTENSION
"""

DEFAULT_RENAME = {'OPTICAL_ELEMENT_NAME':'FILTER', 'EXPOSURE_DURATION':'EXPTIME', 'INSTRUMENT_NAME':'INSTRUMENT', 'DETECTOR_NAME':'DETECTOR', 'STC_S':'FOOTPRINT', 'SPATIAL_RESOLUTION':'PIXSCALE', 'TARGET_NAME':'TARGET', 'SET_ID':'VISIT'}

def demo():
    
    import esa_query
    
    # Query around GN-z11 program
    tab = esa_query.run_query(box=None, proposid=[13871], instruments=['WFC3', 'ACS'], filters=[])
    tab['orientat'] = [esa_query.get_orientat(p) for p in tab['footprint']]
    
    # Search imaging around the center
    box = [np.median(tab['ra']), np.median(tab['dec']), 8]
    
    fig = plt.figure(figsize=[12,3.5])
    ax = fig.add_subplot(141)
    esa_query.show_footprints(tab, ax=ax)
    ax.set_title('13871')
    
    for i, filt in enumerate(['G102', 'F140W', 'F775W']):
        print(filt)
        itab = esa_query.run_query(box=box, proposid=[], instruments=['WFC3', 'ACS'], filters=[filt])
    
        ax = fig.add_subplot(142+i)
        esa_query.show_footprints(itab, ax=ax)
        ax.set_title(filt)
    
    for ax in fig.axes:        
        ax.set_xlim(189.46939409999999, 188.89701790000001)
        ax.set_ylim(62.078094353454389, 62.402351754992722)
        ax.grid()
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    
    fig.tight_layout(pad=0.2)
    fig.savefig('esa_query_demo.png')
    
    # Bielby
    tab = esa_query.run_query(box=None, proposid=[14594], instruments=['WFC3', 'ACS'], extensions=['FLT'], filters=['G141'])
    
    tab = esa_query.run_query(box=None, proposid=[12177], instruments=['WFC3', 'ACS'], extensions=['FLT'], filters=['G141'])
    targets = np.unique(tab['target'])
    
    # WISP
    tab = esa_query.run_query(box=None, proposid=[14178], instruments=['WFC3', 'ACS'], extensions=['FLT'], filters=['G102','G141'])
    targets = np.unique(tab['jtargname'])
    
    for target in targets:
        print(target)
        m = tab['jtargname'] == target
        box = [np.median(tab['ra'][m]), np.median(tab['dec'][m]), 10]
        xtab = esa_query.run_query(box=box, proposid=[], instruments=['WFC3', 'ACS'], extensions=['FLT'], filters=[])
        
        fig = plt.figure(figsize=[5,5])
        ax = fig.add_subplot(111)
        esa_query.show_footprints(tab[m], ax=ax)
        esa_query.show_footprints(xtab, ax=ax)
        ax.set_title(target)
        ax.grid()
        fig.tight_layout(pad=0.1)
        
        fig.savefig('{0}_footprint.png'.format(target))
        
    #
    from shapely.geometry import Polygon
    polygons = []
    for i in range(len(tab)):
        poly = parse_polygons(tab['footprint'][i])#[0]
        pshape = [Polygon(p) for p in poly]
        for i in range(1,len(poly)):
            pshape[0] = pshape[0].union(pshape[i])
        
        polygons.append(pshape[0].buffer(2./60))
    
    match_poly = [polygons[0]]
    match_ids = [[0]]
    
    for i in range(1,len(tab)):
        print(i)
        has_match = False
        for j in range(len(match_poly)):
            isect = match_poly[j].intersection(polygons[i])
            #print(tab['target'][i], i, isect.area > 0)
            if isect.area > 0:
                match_poly[j] = match_poly[j].union(polygons[i])
                match_ids[j].append(i)
                has_match = True
                continue
                
        if not has_match:
            match_poly.append(polygons[i])
            match_ids.append([i])
    
    from descartes import PolygonPatch
    BLUE = '#6699cc'
    for i in range(len(match_poly)):
        #i+=1
        p = match_poly[i]
        print(i)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        idx = np.array(match_ids[i])
        esa_query.show_footprints(tab[idx], ax=ax)

        if True:
            box = [np.median(tab['ra'][idx]), np.median(tab['dec'][idx]), 6]
            xtab = esa_query.run_query(box=box, proposid=[], instruments=['WFC3', 'ACS', 'WFPC2'], extensions=['FLT','C1M'], filters=[])
            esa_query.show_footprints(xtab, ax=ax)
        
        patch1 = PolygonPatch(p, fc=BLUE, ec=BLUE, alpha=0.1, zorder=2)
        
        xy = p.boundary.xy
        ax.plot(xy[0], xy[1])
        ax.add_patch(patch1)
        
        ax.grid()
        
                
def run_query(box=None, proposid=[13871], instruments=['WFC3'], filters=[], extensions=['RAW'], extra=["OBSERVATION.INTENT LIKE 'Science'"],  fields=','.join(DEFAULT_FIELDS.split()), maxitems=1000, rename_columns=DEFAULT_RENAME, lower=True, sort=['OBSERVATION_ID'], remove_tempfile=True):
    """
    
    Optional position box query:
        box = [ra, dec, radius] with ra and dec in decimal degrees and radius
        in arcminutes.
    
    box=None; proposid=[13871]; instruments=['WFC3']; filters=[]; extensions=['RAW']; extra=[];  fields='OBSERVATION,TARGET,POSITION.RA,POSITION.DEC,POSITION.ECL_LAT,POSITION.ECL_LON,POSITION.GAL_LAT,POSITION.GAL_LON,POSITION.STS_S,POSITION.FOV_SIZE,INSTRUMENT.INSTRUMENT_NAME,OPTICAL_ELEMENT.OPTICAL_ELEMENT_NAME,PROPOSAL.PROPOSAL_ID,PROPOSAL.SCIENCE_CATEGORY,PROPOSAL.PI_NAME,ARTIFACT.FILE_NAME'; maxitems=1000
     
    """
    import os
    import tempfile   
    import urllib.request
    from astropy.table import Table
    
    qlist = []
    
    # Box search around position
    if (box is not None):
        ra, dec, radius = box
        dra, ddec = radius/60./np.cos(dec/180*np.pi), radius/60.
        
        bbox = 'POSITION.RA > {0} AND POSITION.RA < {1} AND POSITION.DEC > {2} AND POSITION.DEC < {3}'.format(ra-dra, ra+dra, dec-ddec, dec+ddec)
        
        qlist.append(bbox)
    
    if len(proposid) > 0:
        pquery = ' OR '.join(['PROPOSAL.PROPOSAL_ID=={0}'.format(p) for p in proposid])
        qlist.append('({0})'.format(pquery))
    
    if len(instruments) > 0:
        iquery = ' OR '.join(['INSTRUMENT.INSTRUMENT_NAME LIKE \'{0}\''.format(p) for p in instruments])
        qlist.append('({0})'.format(iquery))
    
    if len(filters) > 0:
        fquery = ' OR '.join(['OPTICAL_ELEMENT.OPTICAL_ELEMENT_NAME LIKE \'{0}\''.format(p) for p in filters])
        qlist.append('({0})'.format(fquery))
    
    if len(extensions) > 0:
        equery = ' OR '.join(['ARTIFACT.FILE_EXTENSION LIKE \'{0}\''.format(p) for p in extensions])
        qlist.append('({0})'.format(equery))
    
            
    query = "http://archives.esac.esa.int/ehst-sl-server/servlet/metadata-action?RESOURCE_CLASS=OBSERVATION&QUERY=({0})&SELECTED_FIELDS={1}&PAGE=1&PAGE_SIZE={2}&RETURN_TYPE=CSV".format(' AND '.join(qlist+extra), fields, maxitems).replace(' ','%20')
    
    req = urllib.request.Request(query)
    response = urllib.request.urlopen(req)
    the_page = response.read().decode('utf-8')
    
    if len(the_page) == 0:
        print('Empty query: ')
        print('\n', '\n '.join(qlist))
        print('\n', fields)
        return False
        
    fp = tempfile.NamedTemporaryFile('w', delete=False)
    fp.write(the_page.replace('"',''))
    fp.close()

    tab = Table.read(fp.name, format='csv')

    if remove_tempfile:
        os.unlink(fp.name)
    else:
        print('Temporary CSV file: ', fp.name)

    # 
    # 
    # header = the_page.split('\n')[0]
    # cols = header.split(',')  
    #tab = Table.read(the_page.split('\n'), format='csv')
        
    # Sort by observation ID
    
    tab.sort(sort)
    
    # Add coordinate name
    if 'RA' in tab.colnames:
        jtargname = [radec_to_targname(ra=tab['RA'][i], dec=tab['DEC'][i], scl=6) for i in range(len(tab))]
        tab['JTARGNAME'] = jtargname
        
    for c in rename_columns:
        if c in tab.colnames:
            tab.rename_column(c, rename_columns[c])
    
    if lower:
        for c in tab.colnames:
            tab.rename_column(c, c.lower())
            
    #tab['OBSERVATION_ID','orientat'][so].show_in_browser(jsviewer=True)
    
    return tab

def radec_to_targname(ra=0, dec=0, scl=10000):
    """Turn decimal degree coordinates into a string
    
    Example:

        >>> from grizli.utils import radec_to_targname
        >>> print(radec_to_targname(ra=10., dec=-10.))
        j004000-100000
    
    Parameters
    -----------
    ra, dec : float
        Sky coordinates in decimal degrees
    
    Returns
    --------
    targname : str
        Target name like jHHMMSS[+-]DDMMSS.
    
    """
    import astropy.coordinates 
    import astropy.units as u
    import re
    import numpy as np
    
    dec_scl = int(np.round(dec*scl))/scl
    cosd = np.cos(dec_scl/180*np.pi)
    ra_scl = int(np.round(ra*scl/cosd))/(scl/cosd)
    
    coo = astropy.coordinates.SkyCoord(ra=ra_scl*u.deg, dec=dec_scl*u.deg)
    
    cstr = re.split('[hmsd.]', coo.to_string('hmsdms', precision=2))
    targname = ('j{0}{1}'.format(''.join(cstr[0:3]), ''.join(cstr[4:7])))
    targname = targname.replace(' ', '')
    
    return targname
        
def parse_polygons(polystr):
    if 'UNION' in polystr.upper():
        spl = polystr[:-1].split('Polygon')[1:]
    else:
        spl = polystr.replace('FK5','ICRS').split('ICRS')[1:]
        
    poly = [np.cast[float](p.split()).reshape((-1,2)) for p in spl]
    return poly
    
def get_orientat(polystr='Polygon ICRS 127.465487 18.855605 127.425760 18.853486 127.423118 18.887458 127.463833 18.889591'):
    """
    
    Compute the "ORIENTAT" position angle (PA of the detector +y axis) from an 
    ESA archive polygon
    
    """
    from astropy.coordinates import Angle
    import astropy.units as u
    
    p = parse_polygons(polystr)[0]
    
    dra = (p[1,0]-p[0,0])*np.cos(p[0,1]/180*np.pi)
    dde = p[1,1] - p[0,1]
    
    orientat = 90+np.arctan2(dra, dde)/np.pi*180
    orientat -= 0.24 # small offset to better match header keywords
    
    orientat = Angle.wrap_at(orientat*u.deg, 180*u.deg).value
    
    return orientat
    
def show_footprints(tab, ax=None):
    import matplotlib.pyplot as plt
    
    # Show polygons
    mpl_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    filters = np.unique(tab['filter'])
    #print(filters)
    colors = {}
    
    if ax is None:
        ax = plt
        
    for i, f in enumerate(filters):
        #print(f)
        if f in MASTER_COLORS:
            colors[f] = MASTER_COLORS[f]
        else:
            colors[f] = mpl_colors[i % len(mpl_colors)]
        
    for i in range(len(tab)):
        poly = parse_polygons(tab['footprint'][i])#[0]
        #print('x',poly, tab['footprint'][i])
        for p in poly:
            #print(p, p.shape)
            pclose = np.vstack([p, p[0,:]])
            #print(pclose)
            ax.plot(pclose[:,0], pclose[:,1], alpha=0.1, color=colors[tab['filter'][i]])
            ax.scatter(pclose[0,0], pclose[0,1], marker='.', color=colors[tab['filter'][i]], alpha=0.1)
    
    #plt.plot(p[0::2], p[1::2], alpha=0.5, color='r')
    #plt.scatter(p[0], p[1], marker='o', color='r')
    
    # orientat = PA y, first two elements of polygon are (y @ x=0)

