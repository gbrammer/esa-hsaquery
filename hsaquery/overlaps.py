"""
Scripts to find overlapping HST data
"""

try:
    from . import query, utils
    from .query import parse_polygons
except:
    from hsaquery import query, utils
    from hsaquery.query import parse_polygons
    
def test():
    
    import copy
    
    import numpy as np
    import matplotlib.pyplot as plt

    from shapely.geometry import Polygon
    from descartes import PolygonPatch
    
    from hsaquery import query
    from hsaquery.query import parse_polygons

    
    # Example: high-z cluster pointings
    tab = query.run_query(box=None, proposid=[14594], instruments=['WFC3', 'ACS'], extensions=['FLT'], filters=['G102','G141'], extra=[])
    
    # Ebeling
    tab = query.run_query(box=None, proposid=[15132,14098,13671,12884,12166,11103,10875,], instruments=['WFC3'], extensions=['FLT'], filters=[], extra=[])
    # Relics
    tab = query.run_query(box=None, proposid=[14096], instruments=['WFC3'], extensions=['FLT'], filters=[], extra=[])
    tab = query.run_query(box=None, proposid=[11591], instruments=['WFC3'], extensions=['FLT'], filters=[], extra=[])
    tab = query.run_query(box=None, proposid=[13666,14148,14496], instruments=['WFC3'], extensions=['FLT'], filters=[], extra=[])
    tab = tab[tab['target'] != 'ANY']
    
    # MACS 0454
    box = [73.5462181, -3.0147200, 3]
    tab = query.run_query(box=box, proposid=[], instruments=['WFC3', 'ACS'], extensions=['FLT'], filters=['F110W'], extra=[])
    
def find_overlaps(tab, buffer_arcmin=1., filters=[], instruments=['WFC3', 'ACS'], proposid=[], SKIP=False, extra=query.DEFAULT_EXTRA, close=True):
    
    import copy
    import os
    
    import numpy as np
    import matplotlib.pyplot as plt

    from shapely.geometry import Polygon
    from descartes import PolygonPatch
        
    # Get shapely polygons for each exposures
    polygons = []
    
    poly_buffer = buffer_arcmin/60 # ~1 arcmin, but doesn't account for cos(dec)
    #poly_buffer = 0.5/60 # ~1 arcmin, but doesn't account for cos(dec)
    
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
        print('Parse', i)
        has_match = False
        for j in range(len(match_poly)):
            isect = match_poly[j].intersection(polygons[i])
            #print(tab['target'][i], i, isect.area > 0)
            if isect.area*3600. > 0.5:
                #print(isect.area*3600)
                match_poly[j] = match_poly[j].union(polygons[i])
                match_ids[j].append(i)
                has_match = True
                continue
                
        if not has_match:
            match_poly.append(polygons[i])
            match_ids.append([i])
    
    ##################
    # Iterate joining polygons
    for iter in range(3):
        mpolygons = copy.deepcopy(match_poly)
        mids = copy.deepcopy(match_ids)
    
        match_poly = [mpolygons[0]]
        match_ids = [mids[0]]
    
        for i in range(1,len(mpolygons)):
            #print(i)
            has_match = False
            for j in range(len(match_poly)):
                isect = match_poly[j].intersection(mpolygons[i])
                #print(tab['target'][i], i, isect.area > 0)
                if isect.area > 0:
                    print(isect.area*3600)
                    match_poly[j] = match_poly[j].union(mpolygons[i])
                    match_ids[j].extend(mids[i])
                    has_match = True
                    continue

            if not has_match:
                match_poly.append(mpolygons[i])
                match_ids.append(mids[i])
        
        print('Iter #{0}, N_Patch = {1}'.format(iter+1, len(match_poly)))
        if len(mpolygons) == len(match_poly):
            break
            
    # Save figures and tables for the unique positions
    BLUE = '#6699cc'
    
    tables = []
    
    for i in range(len(match_poly)):
        #i+=1
        p = match_poly[i]
                
        #######
        # Query around central coordinate
        idx = np.array(match_ids[i])
        ra, dec = np.mean(tab['ra'][idx]), np.mean(tab['dec'][idx])
        
        # Get poly size
        xy = p.convex_hull.boundary.xy
        xradius = np.abs(xy[0]-ra).max()*np.cos(dec/180*np.pi)*60
        yradius = np.abs(xy[1]-dec).max()*60

        box = [ra, dec, np.maximum(xradius*1.5, yradius*1.5)]
        
        # Build target name from RA/Dec
        jname = utils.radec_to_targname(box[0], box[1], scl=1000)
        print('\n\n', i, jname, box[0], box[1])

        if (os.path.exists('{0}_footprint.pdf'.format(jname))) & SKIP:
            continue
                            
        xtab = query.run_query(box=box, proposid=proposid, instruments=instruments, extensions=['FLT','C1M'], filters=filters, extra=extra)
        
        ebv = utils.get_irsa_dust(ra, dec, type='SandF')
        xtab.meta['NAME'] = jname
        xtab.meta['RA'] = ra
        xtab.meta['DEC'] = dec
        xtab.meta['MW_EBV'] = ebv
              
        # Only include ancillary data that directly overlaps with the primary
        # polygon
        pointing_overlaps = np.zeros(len(xtab), dtype=bool)
        for j in range(len(xtab)):
            poly = parse_polygons(xtab['footprint'][j])#[0]
            pshape = [Polygon(pi) for pi in poly]
            for k in range(1,len(poly)):
                pshape[0] = pshape[0].union(pshape[k])
            
            isect = p.intersection(pshape[0])
            pointing_overlaps[j] = isect.area > 0
            
        xtab = xtab[pointing_overlaps]
        
        # Unique targets
        filter_target = np.array(['{0} {1}'.format(xtab['instdet'][i], xtab['filter'][i]) for i in range(pointing_overlaps.sum())])
                        
        ########### 
        # Make the figure
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        # Show the parent table
        colors = query.show_footprints(tab[idx], ax=ax)
        ax.scatter(box[0], box[1], marker='+', color='k')
        
        colors = query.show_footprints(xtab, ax=ax)
        
        patch1 = PolygonPatch(p, fc=BLUE, ec=BLUE, alpha=0.1, zorder=2)
        
        ax.plot(xy[0], xy[1])
        ax.add_patch(patch1)
        
        ax.grid()
        
        # Resize for square dimensions
        xr, yr = ax.get_xlim(), ax.get_ylim()
        dx = (xr[1]-xr[0])*np.cos(yr[0]/180*np.pi)*60
        dy = (yr[1]-yr[0])*60
        ax.set_title(jname)
        ax.set_xlim(ax.get_xlim()[::-1])
        fig.set_size_inches(5,5*dy/dx)
        
        # Add summary
        fp = open('{0}_info.dat'.format(jname),'w')
        dyi = 0.02*dx/dy
        
        ax.text(0.05, 0.97, '{0:>13.5f} {1:>13.5f}  E(B-V)={2:.3f}'.format(ra, dec, ebv), ha='left', va='top', transform=ax.transAxes, fontsize=6)
        
        for i, t in enumerate(np.unique(xtab['proposal_id'])):
            fp.write('proposal_id {0} {1}\n'.format(jname, t))
            ts = np.unique(xtab['target'][xtab['proposal_id'] == t])
            ax.text(0.05, 0.97-dyi*(i+1), '{0} {1}'.format(t, ' '.join(['{0}'.format(ti) for ti in ts])), ha='left', va='top', transform=ax.transAxes, fontsize=6)
            
        for i, t in enumerate(np.unique(xtab['target'])):
            fp.write('target {0} {1}\n'.format(jname, t))
            
        print(np.unique(xtab['target']), '\n')
        
        for i, filt in enumerate(np.unique(filter_target)):
            mf = filter_target == filt
            print('filter {0}  {1:>20s}  {2:>3d}  {3:>8.1f}'.format(jname, filt, mf.sum(), xtab['exptime'][mf].sum()))
            fp.write('filter {0}  {1:>20s}  {2:>3d}  {3:>8.1f}\n'.format(jname, filt, mf.sum(), xtab['exptime'][mf].sum()))
            
            c = colors[filt.split()[1]]
            ax.text(0.95, 0.97-dyi*i, '{1:>20s}  {2:>3d}  {3:>8.1f}\n'.format(jname, filt, mf.sum(), xtab['exptime'][mf].sum()), ha='right', va='top', transform=ax.transAxes, fontsize=6, color=c)
            
        fp.close()
                
        fig.tight_layout(pad=0.5)
        
        # Save figure and table
        fig.savefig('{0}_footprint.pdf'.format(jname))
        
        if close:
            plt.close()

        xtab.write('{0}_footprint.fits'.format(jname), format='fits', overwrite=True)
        np.save('{0}_footprint.npy'.format(jname), [p, box])
        
        tables.append(xtab)
    
    return tables
    
    