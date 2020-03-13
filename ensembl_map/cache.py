"""
Select the Ensembl release to use and download/index the annotations if needed.

Todos:
    * User specify release/species
    * User specify cache directory
"""
from pyensembl import EnsemblRelease


DEFAULT_RELEASE = 69
DEFAULT_SPECIES = "homo_sapiens"

data = EnsemblRelease(DEFAULT_RELEASE, species=DEFAULT_SPECIES)

# download missing annotation files, won't download if they already exist
data.download()
data.index()
