"""
N.B. Suppress all warnings in all of kineticsTools.

The main problem is H5pyDeprecationWarning.
To address that, I think we would need to drop ".value" for h5py
calls in pbcore, which I'd rather not touch.

https://jira.pacificbiosciences.com/browse/TAG-4336
"""
import warnings
warnings.filterwarnings("ignore")
