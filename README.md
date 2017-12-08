# esa-hsaquery
Python tools for querying the ESA Hubble Science Archive (http://archives.esac.esa.int/ehst/#search)

## Installation:

    pip install hsaquery
    
## Demo:

```python
>>> from hsaquery import query

>>> tab = query.run_query(box=None, proposid=[11359], instruments=['WFC3'], 
                     extensions=['FLT'], filters=['G141'], extra=[])
                     
>>> print(tab['observation_id', 'filter', 'exptime'])
observation_id  filter  exptime
--------------  ------  -------
     IB6O23RSQ    G141     1103
     IB6O23RUQ    G141     1003
     IB6O23RYQ    G141     1103
     IB6O23S0Q    G141     1003
``` 
