# esa-hsaquery
Python tools for querying the ESA Hubble Science Archive (http://archives.esac.esa.int/ehst/#search)

## Installation:

    # From PIP
    pip install hsaquery
    
    # Latest version of the respository
    git clone https://github.com/gbrammer/esa-hsaquery.git
    cd esa-hsaquery
    python setup.py install
    
## Demo:

```python
>>> from hsaquery import query

>>> tab = query.run_query(box=None, proposid=[11359], instruments=['WFC3-IR'], 
                     extensions=['FLT'], filters=['G141'], extra=query.DEFAULT_EXTRA)

# NB: as of May 2018 the exptime values from the archive
#     are incorrect in that they are the total exposure time 
#     for the visit that contains a given exposure.                     
>>> print(tab['observation_id', 'filter', 'exptime'])
observation_id filter exptime
-------------- ------ -------
     ib6o23rsq   G141    4212
     ib6o23ruq   G141    4212
     ib6o23ryq   G141    4212
     ib6o23s0q   G141    4212
``` 
