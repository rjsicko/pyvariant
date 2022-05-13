from __future__ import annotations

from dataclasses import dataclass
from locale import normalize
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import numpy
import pandas
from logzero import logger

from ensembl_map.data_cache.data_cache import DataCache
from ensembl_map.utils import is_ensembl_id

from .constants import CDS, CONTIG, EXON, GENE, PROTEIN, TRANSCRIPT

# from .ensembl import Genome, load_ensembl
# from .utils import strip_version

DEFAULT_SPECIES = "homo_sapiens"
DEFAULT_REFERENCE = "GRCh38"
DEFAULT_RELEASE = 100

# column names
# ['gene_id', 'gene_version', 'gene_name', 'gene_source', 'gene_biotype',
# 'transcript_id', 'transcript_version', 'transcript_name', 'transcript_source',
# 'transcript_biotype', 'tag', 'transcript_support_level', 'exon_number', 'exon_id',
# 'exon_version', 'protein_id', 'protein_version', 'ccds_id']
DEFAULT_COLS = [
    "seqname",
    "start",
    "end",
    "strand",
]
EXON_INFO_COLS = DEFAULT_COLS + [
    "gene_id",
    "gene_version",
    "gene_name",
    "gene_biotype",
    "transcript_id",
    "transcript_version",
    "transcript_name",
    "transcript_biotype",
    "exon_number",
    "exon_id",
    "exon_version",
    "protein_id",
    "protein_version",
]
GENE_INFO_COLS = DEFAULT_COLS + [
    "gene_id",
    "gene_version",
    "gene_name",
    "gene_biotype",
]
TRANSCRIPT_INFO_COLS = DEFAULT_COLS + [
    "gene_id",
    "gene_version",
    "gene_name",
    "gene_biotype",
    "transcript_id",
    "transcript_version",
    "transcript_name",
    "transcript_biotype",
    "protein_id",
    "protein_version",
]


class Record:
    def __init__(self, items: Dict):
        self.__dict__.update(items)


class CachedEnsemblRelease(type):
    _instances: Dict[Any, EnsemblRelease] = {}
    _current: Optional[EnsemblRelease] = None

    def __call__(cls, *args, **kwargs):
        # convert the given arguments to a hashable tuple
        key = tuple(list(args) + [kwargs[k] for k in sorted(kwargs)])
        # check if an instance created with the same arguments is already cached
        if key in cls._instances:
            instance = cls._instances[key]
        else:
            # create a new instance and cache it
            instance = super(CachedEnsemblRelease, cls).__call__(*args, **kwargs)
            cls._instances[key] = instance

        # track which instance was called last
        cls._current = instance

        return cls._instances[key]


class EnsemblRelease(metaclass=CachedEnsemblRelease):
    def __init__(
        self,
        species: str = DEFAULT_SPECIES,
        reference: str = DEFAULT_REFERENCE,
        release: int = DEFAULT_RELEASE,
        cache_dir: str = "",
        canonical_transcript: str = "",
        contig_alias: str = "",
        gene_alias: str = "",
        transcript_alias: str = "",
    ):
        """Load annotations for the given release."""
        self.species = species
        self.reference = reference
        self.release = release
        self.cache_dir = cache_dir

        self.data_cache = DataCache(species, reference, release, cache_dir=cache_dir).load_gtf()

        self.canonical_transcript = _parse_txt_to_list(canonical_transcript, "canonical transcript")
        self.contig_alias = _parse_tsv_to_dict(contig_alias, "contig aliases")
        self.exon_alias: Dict[str, List[str]] = {}
        self.gene_alias = _parse_tsv_to_dict(gene_alias, "gene aliases")
        self.protein_alias: Dict[str, List[str]] = {}
        self.transcript_alias = _parse_tsv_to_dict(transcript_alias, "transcript aliases")

    # ---------------------------------------------------------------------------------------------
    # <feature_symbol>()
    # ---------------------------------------------------------------------------------------------
    def all_contig_ids(self) -> List[str]:
        return self._uniquify_series(self.data_cache["seqname"])

    def all_exon_ids(self) -> List[str]:
        return self._uniquify_series(self.data_cache["exon_id"])

    def all_gene_ids(self) -> List[str]:
        return self._uniquify_series(self.data_cache["gene_id"])

    def all_gene_names(self) -> List[str]:
        return self._uniquify_series(self.data_cache["gene_name"])

    def all_protein_ids(self) -> List[str]:
        return self._uniquify_series(self.data_cache["protein_id"])

    def all_transcript_ids(self) -> List[str]:
        return self._uniquify_series(self.data_cache["transcript_id"])

    def all_transcript_names(self) -> List[str]:
        return self._uniquify_series(self.data_cache["transcript_name"])

    # ---------------------------------------------------------------------------------------------
    # get_<feature_symbol>(feature)
    # ---------------------------------------------------------------------------------------------
    def contig_ids(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding contig ID(s)."""
        result = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._contig_ids_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._contig_ids_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._contig_ids_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._contig_ids_of_protein_id(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._contig_ids_of_transcript_id(feature))
                else:
                    raise ValueError(f"Unable to get contig IDs for {feature} ({feature_type})")
            else:
                if feature_type == CONTIG:
                    result.extend(self._contig_ids_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._contig_ids_of_gene_name(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._contig_ids_of_transcript_name(feature))
                else:
                    raise ValueError(f"Unable to get contig IDs for {feature} ({feature_type})")

        return result

    def exon_ids(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding exon ID(s)."""
        result = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._exon_ids_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._exon_ids_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._exon_ids_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._exon_ids_of_protein_id(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._exon_ids_of_transcript_id(feature))
                else:
                    raise ValueError(f"Unable to get exon IDs for {feature} ({feature_type})")
            else:
                if feature_type == CONTIG:
                    result.extend(self._exon_ids_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._exon_ids_of_gene_name(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._exon_ids_of_transcript_name(feature))
                else:
                    raise ValueError(f"Unable to get exon IDs for {feature} ({feature_type})")

        return result

    def gene_ids(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding gene ID(s)."""
        result = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._gene_ids_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._gene_ids_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._gene_ids_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._gene_ids_of_protein_id(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._gene_ids_of_transcript_id(feature))
                else:
                    raise ValueError(f"Unable to get gene IDs for {feature} ({feature_type})")
            else:
                if feature_type == CONTIG:
                    result.extend(self._gene_ids_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._gene_ids_of_gene_name(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._gene_ids_of_transcript_name(feature))
                else:
                    raise ValueError(f"Unable to get gene IDs for {feature} ({feature_type})")

        return result

    def gene_names(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding gene names(s)."""
        result = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._gene_names_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._gene_names_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._gene_names_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._gene_names_of_protein_id(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._gene_names_of_transcript_id(feature))
                else:
                    raise ValueError(f"Unable to get gene names for {feature} ({feature_type})")
            else:
                if feature_type == CONTIG:
                    result.extend(self._gene_names_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._gene_names_of_gene_name(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._gene_names_of_transcript_name(feature))
                else:
                    raise ValueError(f"Unable to get gene names for {feature} ({feature_type})")

        return result

    def protein_ids(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding protein ID(s)."""
        result = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._protein_ids_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._protein_ids_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._protein_ids_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._protein_ids_of_protein_id(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._protein_ids_of_transcript_id(feature))
                else:
                    raise ValueError(f"Unable to get protein IDs for {feature} ({feature_type})")
            else:
                if feature_type == CONTIG:
                    result.extend(self._protein_ids_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._protein_ids_of_gene_name(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._protein_ids_of_transcript_name(feature))
                else:
                    raise ValueError(f"Unable to get protein IDs for {feature} ({feature_type})")

        return result

    def transcript_ids(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding transcript ID(s)."""
        result = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._transcript_ids_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._transcript_ids_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._transcript_ids_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._transcript_ids_of_protein_id(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._transcript_ids_of_transcript_id(feature))
                else:
                    raise ValueError(f"Unable to get transcript IDs for {feature} ({feature_type})")
            else:
                if feature_type == CONTIG:
                    result.extend(self._transcript_ids_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._transcript_ids_of_gene_name(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._transcript_ids_of_transcript_name(feature))
                else:
                    raise ValueError(f"Unable to get transcript IDs for {feature} ({feature_type})")

        return result

    def transcript_names(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding transcript names(s)."""
        result = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._transcript_names_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._transcript_names_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._transcript_names_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._transcript_names_of_protein_id(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._transcript_names_of_transcript_id(feature))
                else:
                    raise ValueError(
                        f"Unable to get transcript names for {feature} ({feature_type})"
                    )
            else:
                if feature_type == CONTIG:
                    result.extend(self._transcript_ids_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._transcript_names_of_gene_name(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._transcript_names_of_transcript_name(feature))
                else:
                    raise ValueError(
                        f"Unable to get transcript names for {feature} ({feature_type})"
                    )

        return result

    def _contig_ids_of_contig_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "seqname", "seqname")

    def _contig_ids_of_exon_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "exon_id", "seqname")

    def _contig_ids_of_gene_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_id", "seqname")

    def _contig_ids_of_gene_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_name", "seqname")

    def _contig_ids_of_protein_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "protein_id", "seqname")

    def _contig_ids_of_transcript_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_id", "seqname")

    def _contig_ids_of_transcript_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_name", "seqname")

    def _exon_ids_of_contig_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "seqname", "exon_id")

    def _exon_ids_of_exon_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "exon_id", "exon_id")

    def _exon_ids_of_gene_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_id", "exon_id")

    def _exon_ids_of_gene_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_name", "exon_id")

    def _exon_ids_of_protein_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(self._transcript_ids_of_protein_id(feature), "transcript_id", "exon_id")

    def _exon_ids_of_transcript_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_id", "exon_id")

    def _exon_ids_of_transcript_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_name", "exon_id")

    def _gene_ids_of_contig_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "seqname", "gene_id")

    def _gene_ids_of_exon_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "exon_id", "gene_id")

    def _gene_ids_of_gene_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_id", "gene_id")

    def _gene_ids_of_gene_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_name", "gene_id")

    def _gene_ids_of_protein_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "protein_id", "gene_id")

    def _gene_ids_of_transcript_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_id", "gene_id")

    def _gene_ids_of_transcript_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_name", "gene_id")

    def _gene_names_of_contig_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "seqname", "gene_name")

    def _gene_names_of_exon_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "exon_id", "gene_name")

    def _gene_names_of_gene_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_id", "gene_name")

    def _gene_names_of_gene_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_name", "gene_name")

    def _gene_names_of_protein_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "protein_id", "gene_name")

    def _gene_names_of_transcript_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_id", "gene_name")

    def _gene_names_of_transcript_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_name", "gene_name")

    def _protein_ids_of_contig_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "seqname", "protein_id")

    def _protein_ids_of_exon_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(self._transcript_ids_of_exon_id(feature), "transcript_id", "protein_id")

    def _protein_ids_of_gene_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_id", "protein_id")

    def _protein_ids_of_gene_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_name", "protein_id")

    def _protein_ids_of_protein_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "protein_id", "protein_id")

    def _protein_ids_of_transcript_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_id", "protein_id")

    def _protein_ids_of_transcript_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_name", "protein_id")

    def _transcript_ids_of_contig_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "seqname", "transcript_id")

    def _transcript_ids_of_exon_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "exon_id", "transcript_id")

    def _transcript_ids_of_gene_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_id", "transcript_id")

    def _transcript_ids_of_gene_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_name", "transcript_id")

    def _transcript_ids_of_protein_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "protein_id", "transcript_id")

    def _transcript_ids_of_transcript_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_id", "transcript_id")

    def _transcript_ids_of_transcript_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_name", "transcript_id")

    def _transcript_names_of_contig_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "seqname", "transcript_name")

    def _transcript_names_of_exon_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "exon_id", "transcript_name")

    def _transcript_names_of_gene_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_id", "transcript_name")

    def _transcript_names_of_gene_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_name", "transcript_name")

    def _transcript_names_of_protein_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "protein_id", "transcript_name")

    def _transcript_names_of_transcript_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_id", "transcript_name")

    def _transcript_names_of_transcript_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_name", "transcript_name")

    def _query(self, feature: Union[List[str], str], col: str, key: str) -> List[str]:
        """Generic function for querying the data cache."""
        feature = [feature] if isinstance(feature, str) else feature
        return self._uniquify_series(self.data_cache.loc[self.data_cache[col].isin(feature)][key])

    @staticmethod
    def _uniquify_series(values: pandas.Series) -> List:
        """Convert a `pandas.Series` object into a sorted list of unique values."""
        return sorted((i for i in values.unique().tolist() if i))

    # ---------------------------------------------------------------------------------------------
    # normalize_feature(feature, feature_type)
    # ---------------------------------------------------------------------------------------------
    def normalize_feature(self, feature: str, feature_type: str = "") -> List[Tuple[str, str]]:
        """Normalize a feature to the representation used by Ensembl."""
        normalized = []

        # first, check if the feature is found in the database
        if feature_type == CONTIG or not feature_type:
            if normalized := self._normalize_contig_id(feature):
                feature_type = CONTIG

        if feature_type == EXON or not feature_type:
            if normalized := self._normalize_exon_id(feature):
                feature_type = EXON

        if feature_type == GENE or not feature_type:
            if normalized := self._normalize_gene_id(feature):
                feature_type = GENE
            elif normalized := self._normalize_gene_name(feature):
                feature_type = GENE

        if feature_type == PROTEIN or not feature_type:
            if normalized := self._normalize_protein_id(feature):
                feature_type = PROTEIN

        if feature_type in (CDS, TRANSCRIPT) or not feature_type:
            if normalized := self._normalize_transcript_id(feature):
                feature_type = TRANSCRIPT
            elif normalized := self._normalize_transcript_name(feature):
                feature_type = TRANSCRIPT

        return [(i, feature_type) for i in normalized]

    def _normalize_contig_id(self, feature: str) -> List[str]:
        """Normalize a contig ID or it's alias to one or more matching Ensembl contig ID."""
        featurel = [feature] + self.contig_alias.get(feature, [])
        if normalized := self._contig_ids_of_contig_id(featurel):
            return normalized
        else:
            return []

    def _normalize_exon_id(self, feature: str) -> List[str]:
        """Normalize a exon ID or it's alias to one or more matching Ensembl exon ID."""
        featurel = [feature] + self.exon_alias.get(feature, [])
        if normalized := self._exon_ids_of_exon_id(featurel):
            return normalized
        else:
            return []

    def _normalize_gene_id(self, feature: str) -> List[str]:
        """Normalize a gene ID or it's alias to one or more matching Ensembl gene ID."""
        featurel = [feature] + self.gene_alias.get(feature, [])
        if normalized := self._gene_ids_of_gene_id(featurel):
            return normalized
        else:
            return []

    def _normalize_gene_name(self, feature: str) -> List[str]:
        """Normalize a gene ID or it's alias to one or more matching Ensembl gene ID."""
        featurel = [feature] + self.gene_alias.get(feature, [])
        if normalized := self._gene_names_of_gene_name(featurel):
            return normalized
        else:
            return []

    def _normalize_protein_id(self, feature: str) -> List[str]:
        """Normalize a protein ID or it's alias to one or more matching Ensembl protein ID."""
        featurel = [feature] + self.protein_alias.get(feature, [])
        if normalized := self._protein_ids_of_protein_id(featurel):
            return normalized
        else:
            return []

    def _normalize_transcript_id(self, feature: str) -> List[str]:
        """Normalize a transcript ID or it's alias to one or more matching Ensembl transcript ID."""
        featurel = [feature] + self.transcript_alias.get(feature, [])
        if normalized := self._transcript_ids_of_transcript_id(featurel):
            return normalized
        else:
            return []

    def _normalize_transcript_name(self, feature: str) -> List[str]:
        """Normalize a transcript name or it's alias to one or more matching Ensembl transcript ID."""
        featurel = [feature] + self.transcript_alias.get(feature, [])
        if normalized := self._transcript_names_of_transcript_name(featurel):
            return normalized
        else:
            return []

    # ---------------------------------------------------------------------------------------------
    # is_<feature>(feature)
    # ---------------------------------------------------------------------------------------------
    def is_contig(self, feature: str) -> bool:
        """Return True if the given feature is a contig."""
        return any((i[1] == CONTIG for i in self.normalize_feature(feature, CONTIG)))

    def is_exon(self, feature: str) -> bool:
        """Return True if the given feature is an exon."""
        return any((i[1] == EXON for i in self.normalize_feature(feature, EXON)))

    def is_gene(self, feature: str) -> bool:
        """Return True if the given feature is a gene."""
        return any((i[1] == GENE for i in self.normalize_feature(feature, GENE)))

    def is_protein(self, feature: str) -> bool:
        """Return True if the given feature is a protein."""
        return any((i[1] == PROTEIN for i in self.normalize_feature(feature, PROTEIN)))

    def is_transcript(self, feature: str) -> bool:
        """Return True if the given feature is a transcript."""
        return any((i[1] == TRANSCRIPT for i in self.normalize_feature(feature, TRANSCRIPT)))

    # ---------------------------------------------------------------------------------------------
    # get_<feature>_info>(feature, feature_type)
    # ---------------------------------------------------------------------------------------------
    def contig_info(self, feature: str, feature_type: str = "") -> List[Record]:
        """Given a feature symbol, return information on the corresponding contig(s)."""
        result: List[Record] = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._contig_info_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._contig_info_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._contig_info_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._contig_info_of_protein_id(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._contig_info_of_transcript_id(feature))
                else:
                    raise ValueError(f"Unable to get contig info for {feature} ({feature_type})")
            else:
                if feature_type == CONTIG:
                    result.extend(self._contig_info_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._contig_info_of_gene_name(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._contig_info_of_transcript_name(feature))
                else:
                    raise ValueError(f"Unable to get contig info for {feature} ({feature_type})")

        return result

    def exon_info(self, feature: str, feature_type: str = "") -> List[Record]:
        """Given a feature symbol, return information on the corresponding exon(s)."""
        result = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._exon_info_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._exon_info_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._exon_info_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._exon_info_of_protein_id(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._exon_info_of_transcript_id(feature))
                else:
                    raise ValueError(f"Unable to get exon info for {feature} ({feature_type})")
            else:
                if feature_type == CONTIG:
                    result.extend(self._exon_info_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._exon_info_of_gene_name(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._exon_info_of_transcript_name(feature))
                else:
                    raise ValueError(f"Unable to get exon info for {feature} ({feature_type})")

        return result

    def gene_info(self, feature: str, feature_type: str = "") -> List[Record]:
        """Given a feature symbol, return information on the corresponding gene(s)."""
        result = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._gene_info_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._gene_info_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._gene_info_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._gene_info_of_protein_id(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._gene_info_of_transcript_id(feature))
                else:
                    raise ValueError(f"Unable to get gene info for {feature} ({feature_type})")
            else:
                if feature_type == CONTIG:
                    result.extend(self._gene_info_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._gene_info_of_gene_name(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._gene_info_of_transcript_name(feature))
                else:
                    raise ValueError(f"Unable to get gene info for {feature} ({feature_type})")

        return result

    def protein_info(self, feature: str, feature_type: str = "") -> List[Record]:
        """Given a feature symbol, return information on the corresponding protein(s)."""
        result = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._protein_info_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._protein_info_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._protein_info_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._protein_info_of_protein_id(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._protein_info_of_transcript_id(feature))
                else:
                    raise ValueError(f"Unable to get protein info for {feature} ({feature_type})")
            else:
                if feature_type == CONTIG:
                    result.extend(self._protein_info_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._protein_info_of_gene_name(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._protein_info_of_transcript_name(feature))
                else:
                    raise ValueError(f"Unable to get protein info for {feature} ({feature_type})")

        return result

    def transcript_info(self, feature: str, feature_type: str = "") -> List[Record]:
        """Given a feature symbol, return information on the corresponding transcript(s)."""
        result = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._transcript_info_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._transcript_info_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._transcript_info_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._transcript_info_of_protein_id(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._transcript_info_of_transcript_id(feature))
                else:
                    raise ValueError(
                        f"Unable to get transcript info for {feature} ({feature_type})"
                    )
            else:
                if feature_type == CONTIG:
                    result.extend(self._transcript_info_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._transcript_info_of_gene_name(feature))
                elif feature_type in (CDS, TRANSCRIPT):
                    result.extend(self._transcript_info_of_transcript_name(feature))
                else:
                    raise ValueError(
                        f"Unable to get transcript info for {feature} ({feature_type})"
                    )

        return result

    def _contig_info_of_contig_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_contig_info(feature, "seqname")

    def _contig_info_of_exon_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_contig_info(self._transcript_ids_of_exon_id(feature), "transcript_id")

    def _contig_info_of_gene_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_contig_info(feature, "gene_id")

    def _contig_info_of_gene_name(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_contig_info(feature, "gene_name")

    def _contig_info_of_protein_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_contig_info(self._transcript_ids_of_protein_id(feature), "transcript_id")

    def _contig_info_of_transcript_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_contig_info(feature, "transcript_id")

    def _contig_info_of_transcript_name(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_contig_info(feature, "transcript_name")

    def _exon_info_of_contig_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_exon_info(feature, "seqname")

    def _exon_info_of_exon_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_exon_info(self._transcript_ids_of_exon_id(feature), "transcript_id")

    def _exon_info_of_gene_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_exon_info(feature, "gene_id")

    def _exon_info_of_gene_name(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_exon_info(feature, "gene_name")

    def _exon_info_of_protein_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_exon_info(self._transcript_ids_of_protein_id(feature), "transcript_id")

    def _exon_info_of_transcript_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_exon_info(feature, "transcript_id")

    def _exon_info_of_transcript_name(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_exon_info(feature, "transcript_name")

    def _gene_info_of_contig_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_gene_info(feature, "seqname")

    def _gene_info_of_exon_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_gene_info(self._transcript_ids_of_exon_id(feature), "transcript_id")

    def _gene_info_of_gene_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_gene_info(feature, "gene_id")

    def _gene_info_of_gene_name(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_gene_info(feature, "gene_name")

    def _gene_info_of_protein_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_gene_info(self._transcript_ids_of_protein_id(feature), "transcript_id")

    def _gene_info_of_transcript_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_gene_info(feature, "transcript_id")

    def _gene_info_of_transcript_name(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_gene_info(feature, "transcript_name")

    def _protein_info_of_contig_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_protein_info(feature, "seqname")

    def _protein_info_of_exon_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_protein_info(self._transcript_ids_of_exon_id(feature), "transcript_id")

    def _protein_info_of_gene_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_protein_info(feature, "gene_id")

    def _protein_info_of_gene_name(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_protein_info(feature, "gene_name")

    def _protein_info_of_protein_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_protein_info(
            self._transcript_ids_of_protein_id(feature), "transcript_id"
        )

    def _protein_info_of_transcript_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_protein_info(feature, "transcript_id")

    def _protein_info_of_transcript_name(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_protein_info(feature, "transcript_name")

    def _transcript_info_of_contig_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_transcript_info(feature, "seqname")

    def _transcript_info_of_exon_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_transcript_info(
            self._transcript_ids_of_exon_id(feature), "transcript_id"
        )

    def _transcript_info_of_gene_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_transcript_info(feature, "gene_id")

    def _transcript_info_of_gene_name(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_transcript_info(feature, "gene_name")

    def _transcript_info_of_protein_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_transcript_info(
            self._transcript_ids_of_protein_id(feature), "transcript_id"
        )

    def _transcript_info_of_transcript_id(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_transcript_info(feature, "transcript_id")

    def _transcript_info_of_transcript_name(self, feature: Union[List[str], str]) -> List[Record]:
        return self._query_transcript_info(feature, "transcript_name")

    def _query_contig_info(self, feature: Union[List[str], str], col: str) -> List[Record]:
        """Generic function for querying contig info from the data cache."""
        raise NotImplementedError()  # TODO

    def _query_exon_info(self, feature: Union[List[str], str], col: str) -> List[Record]:
        """Generic function for querying exon info from the data cache."""
        return self._query_info(feature, "exon", col, EXON_INFO_COLS)

    def _query_gene_info(self, feature: Union[List[str], str], col: str) -> List[Record]:
        """Generic function for querying gene info from the data cache."""
        return self._query_info(feature, "transcript", col, GENE_INFO_COLS)

    def _query_protein_info(self, feature: Union[List[str], str], col: str) -> List[Record]:
        """Generic function for querying protein info from the data cache."""
        raise NotImplementedError()  # TODO

    def _query_transcript_info(self, feature: Union[List[str], str], col: str) -> List[Record]:
        """Generic function for querying transcript info from the data cache."""
        return self._query_info(feature, "transcript", col, TRANSCRIPT_INFO_COLS)

    def _query_info(
        self,
        feature: Union[List[str], str],
        feature_col: str,
        id_col: str,
        subset_cols: List[str],
    ) -> List[Record]:
        """Generic function for querying transcript info from the data cache."""
        info = []

        feature = [feature] if isinstance(feature, str) else feature
        select = (self.data_cache["feature"] == feature_col) & (
            self.data_cache[id_col].isin(feature)
        )
        df = self.data_cache.loc[select][subset_cols]
        for _, row in df.iterrows():
            values = Record(row.to_dict())
            if values not in info:
                info.append(values)

        return info


def instance():
    """Return the last called instance of `EnsemblRelease` otherwise initialize an instance with
    the default parameters and return that.
    """
    if EnsemblRelease._current is None:
        EnsemblRelease._current = EnsemblRelease()

    return EnsemblRelease._current


def _parse_txt_to_list(path: str, message: str) -> List[str]:
    """Parse a text file into a list of unique values."""
    result = set()

    if path:
        logger.debug(f"Loading {message} from {path}")
        with open(path, "r") as fh:
            for line in fh:
                if line:
                    result.add(line.strip())
    else:
        logger.debug(f"No {message} to load")

    return sorted(result)


def _parse_tsv_to_dict(path: str, message: str) -> Dict[str, List[str]]:
    """Parse a TSV of one-to-one alias-to-name mappings."""
    result: Dict = {}

    if path:
        # load data from the TSV file, keep only unique values
        logger.debug(f"Loading {message} from {path}")
        with open(path, "r") as fh:
            for line in fh:
                if line:
                    alias, name = line.strip().split("\t")[:2]
                    result.setdefault(alias, set())
                    result[alias].add(name)

        # convert values to sorted lists
        for alias in result:
            result[alias] = sorted(result[alias])
    else:
        logger.debug(f"No {message} to load")

    return result
