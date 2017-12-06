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

DEFAULT_COLUMN_FORMAT = {'start_time_mjd':'.4f',
           'end_time_mjd':'.4f',
           'exptime':'.0f',
           'ra':'.6f',
           'dec':'.6f',
           'ecl_lat':'.6f',
           'ecl_lon':'.6f',
           'gal_lat':'.6f',
           'gal_lon':'.6f',
           'fov_size':'.3f',
           'pixscale':'.3f'}

# Don't get calibrations.  Can't use "INTENT LIKE 'SCIENCE'" because some 
# science observations are flagged as 'Calibration' in the ESA HSA.
DEFAULT_EXTRA = ["TARGET.TARGET_NAME NOT LIKE '{0}'".format(calib) 
                 for calib in ['DARK','EARTH-CALIB', 'TUNGSTEN', 'BIAS',
                               'DARK-EARTH-CALIB', 'DARK-NM', 'DEUTERIUM']]

def run_query(box=None, proposid=[13871], instruments=['WFC3'], filters=[], extensions=['RAW','C1M'], extra=DEFAULT_EXTRA,  fields=','.join(DEFAULT_FIELDS.split()), maxitems=100000, rename_columns=DEFAULT_RENAME, lower=True, sort_column=['OBSERVATION_ID'], remove_tempfile=True):
    """
    
    Optional position box query:
        box = [ra, dec, radius] with ra and dec in decimal degrees and radius
        in arcminutes.
    
    Some science observations are flagged as INTENT = Calibration, so may have
    to run with extra=[] for those cases and strip out true calibs another
    way.
    
    """
    import os
    import tempfile   
    import urllib.request
    from astropy.table import Table
    from . import utils
    
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

    # Sort
    tab.sort(sort_column)
    
    # Add coordinate name
    if 'RA' in tab.colnames:
        jtargname = [utils.radec_to_targname(ra=tab['RA'][i], dec=tab['DEC'][i], scl=6) for i in range(len(tab))]
        tab['JTARGNAME'] = jtargname
        
    for c in rename_columns:
        if c in tab.colnames:
            tab.rename_column(c, rename_columns[c])
    
    if lower:
        for c in tab.colnames:
            tab.rename_column(c, c.lower())
    
    cols = tab.colnames
    if ('instrument' in cols) & ('detector' in cols):
        tab['instdet'] = ['{0}/{1}'.format(tab['instrument'][i], tab['detector'][i]) for i in range(len(tab))]
            
    #tab['OBSERVATION_ID','orientat'][so].show_in_browser(jsviewer=True)
    
    set_default_formats(tab)
    
    return tab

def add_postcard(table, resolution=256):
    
   url = ['http://archives.esac.esa.int/ehst-sl-server/servlet/data-action?OBSERVATION_ID={0}&RETRIEVAL_TYPE=POSTCARD&RESOLUTION={1}'.format(o, resolution) for o in table['observation_id']]
   
   img = ['<a href="{0}"><img src="{0}"></a>'.format(u) for u in url]
   table['postcard'] = img
   
   return True
   
   if False:
       tab = grizli.utils.GTable(table)
       tab['observation_id','filter','orientat','postcard'][tab['visit'] == 1].write_sortable_html('tab.html', replace_braces=True, localhost=True, max_lines=10000, table_id=None, table_class='display compact', css=None)
        
def parse_polygons(polystr):
    if 'UNION' in polystr.upper():
        spl = polystr[:-1].split('Polygon')[1:]
    else:
        spl = polystr.replace('FK5','ICRS').split('ICRS')[1:]
        
    poly = [np.cast[float](p.split()).reshape((-1,2)) for p in spl]
    return poly
           
def set_default_formats(table, formats=DEFAULT_COLUMN_FORMAT):
    """
    Set default print formats
    """
    DEFAULT_FORMATS = {'start_time_mjd':'.4f',
               'end_time_mjd':'.4f',
               'exptime':'.0f',
               'ra':'.6f',
               'dec':'.6f',
               'ecl_lat':'.6f',
               'ecl_lon':'.6f',
               'gal_lat':'.6f',
               'gal_lon':'.6f',
               'fov_size':'.3f',
               'pixscale':'.3f'}
    
    for f in formats:
        if f in table.colnames:
            table[f].format = formats[f]
            
def set_orientat_column(table):
    """
    Make a column in the `table` computing each orientation with
    `get_orientat`.
    """
    table['orientat'] = [query.get_orientat(p) for p in table['footprint']]
    table['orientat'].format = '.1f'
    
def get_orientat(polystr='Polygon ICRS 127.465487 18.855605 127.425760 18.853486 127.423118 18.887458 127.463833 18.889591'):
    """
    
    Compute the "ORIENTAT" position angle (PA of the detector +y axis) from an 
    ESA archive polygon, assuming that the first two entries of the polygon 
    are the LL and UL corners of the detector.
    
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
    """
    Show pointing footprints in a plot
    """
    import matplotlib.pyplot as plt
    
    # Show polygons
    mpl_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    filters = np.unique(tab['filter'])

    colors = {}
    
    if ax is None:
        ax = plt
        
    for i, f in enumerate(filters):
        if f in MASTER_COLORS:
            colors[f] = MASTER_COLORS[f]
        else:
            colors[f] = mpl_colors[i % len(mpl_colors)]
        
    for i in range(len(tab)):
        poly = parse_polygons(tab['footprint'][i])#[0]

        for p in poly:
            pclose = np.vstack([p, p[0,:]]) # repeat the first vertex to close
            
            ax.plot(pclose[:,0], pclose[:,1], alpha=0.1, color=colors[tab['filter'][i]])
            
            # Plot a point at the first vertex, pixel x=y=0.
            ax.scatter(pclose[0,0], pclose[0,1], marker='.', color=colors[tab['filter'][i]], alpha=0.1)
    
    return colors
    

