"""
Scripts to find overlapping HST data
"""
def test():
    
    import numpy as np
    import matplotlib.pyplot as plt

    from shapely.geometry import Polygon
    from descartes import PolygonPatch
    
    from hsaquery import query
    from hsaquery.query import parse_polygons

    
    # Example: high-z cluster pointings
    tab = query.run_query(box=None, proposid=[11597,12203], instruments=['WFC3', 'ACS'], extensions=['FLT'], filters=['G102','G141'], extra=[])
    
    # Get shapely polygons for each exposures
    polygons = []
    
    poly_buffer = 1./60 # ~1 arcmin, but doesn't account for cos(dec)
    for i in range(len(tab)):
        poly = parse_polygons(tab['footprint'][i])#[0]
        pshape = [Polygon(p) for p in poly]
        for i in range(1,len(poly)):
            pshape[0] = pshape[0].union(pshape[i])
        
        polygons.append(pshape[0].buffer(poly_buffer))
    
    match_poly = [polygons[0]]
    match_ids = [[0]]
    
    # Loop through polygons and combine those that overlap
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
    
    # Save figures and tables for the unique positions
    BLUE = '#6699cc'
    
    for i in range(len(match_poly)):
        #i+=1
        p = match_poly[i]
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        idx = np.array(match_ids[i])
        
        # Show the parent table
        query.show_footprints(tab[idx], ax=ax)
        
        #######
        # Query around central coordinate
        box = [np.mean(tab['ra'][idx]), np.mean(tab['dec'][idx]), 6]
                
        xtab = query.run_query(box=box, proposid=[], instruments=['WFC3', 'ACS'], extensions=['FLT','C1M'], filters=[], extra=["TARGET.TARGET_NAME NOT LIKE '{0}'".format(calib) for calib in ['DARK','EARTH-CALIB', 'TUNGSTEN', 'BIAS', 'DARK-EARTH-CALIB', 'DARK-NM', 'DEUTERIUM']])
        
        # Build target name from RA/Dec
        jname = query.radec_to_targname(box[0], box[1], scl=1000)
        print('\n\n', i, jname, box[0], box[1])
        ax.scatter(box[0], box[1], marker='+', color='k')
        
        # Only include ancillary data that directly overlaps with the primary
        # polygon
        pointing_overlaps = np.zeros(len(xtab), dtype=bool)
        for j in range(len(xtab)):
            poly = parse_polygons(xtab['footprint'][j])#[0]
            pshape = [Polygon(p) for p in poly]
            for k in range(1,len(poly)):
                pshape[0] = pshape[0].union(pshape[k])
            
            isect = p.intersection(pshape[0])
            pointing_overlaps[j] = isect.area > 0
            
        xtab = xtab[pointing_overlaps]
        
        # Unique targets
        filter_target = np.array(['{0} {1}'.format(xtab['instdet'][i], xtab['filter'][i]) for i in range(pointing_overlaps.sum())])
        print(np.unique(xtab['target']), '\n')
        for filt in np.unique(filter_target):
            mf = filter_target == filt
            print('    {0:20s}  {1:3d}  {2:>8.1f}'.format(filt, mf.sum(), xtab['exptime'][mf].sum()))
        
        # Add them to the plot
        query.show_footprints(xtab, ax=ax)
        
        patch1 = PolygonPatch(p, fc=BLUE, ec=BLUE, alpha=0.1, zorder=2)
        
        xy = p.boundary.xy
        ax.plot(xy[0], xy[1])
        ax.add_patch(patch1)
        
        ax.grid()
        
        # Resize for square dimensions
        xr, yr = ax.get_xlim(), ax.get_ylim()
        dx = (xr[1]-xr[0])*np.cos(yr[0]/180*np.pi)*60
        dy = (yr[1]-yr[0])*60
        ax.set_title(jname)
        fig.set_size_inches(5,5*dy/dx)
        fig.tight_layout(pad=0.5)
        
        # Save figure and table
        fig.savefig('{0}_footprint.pdf'.format(jname))
        plt.close()

        xtab.write('{0}_footprint.fits'.format(jname), format='fits', overwrite=True)
        
    