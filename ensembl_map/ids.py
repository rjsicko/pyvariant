from functools import lru_cache
from typing import Any, Callable, List, Optional, Tuple

import pyensembl

from .best_transcripts import BestTranscript
from .constants import CDS, CONTIG, EXON, GENE, PROTEIN, TRANSCRIPT
from .core import instance as CM
from .utils import is_ensembl_id


# --------------------------------------
#  Normalization functions
# --------------------------------------
@lru_cache()
def normalize_feature(feature: str, feature_type: Optional[str] = None) -> List[Tuple[str, str]]:
    normalized = []

    # first, check if the feature is found in the database
    if not feature_type:
        if (
            feature in CM().ensembl.all_transcript_ids()
            or feature in CM().ensembl.all_transcript_names()
        ):
            feature_type = TRANSCRIPT
        elif feature in CM().ensembl.all_gene_ids() or feature in CM().ensembl.all_gene_names():
            feature_type = GENE
        elif feature in CM().ensembl.contigs():
            feature_type = CONTIG
        elif feature in CM().ensembl.all_protein_ids():
            feature_type = PROTEIN
        elif feature in CM().ensembl.all_exon_ids():
            feature_type = EXON

    # second, check if the feature is known by an alias
    if feature_type == CONTIG or not feature_type:
        for contig_id in CM().normalize_contig_id(feature):
            normalized.append((contig_id, CONTIG))

    if feature_type == GENE or not feature_type:
        for gene_id in CM().normalize_gene_id(feature):
            normalized.append((gene_id, GENE))

    if feature_type == TRANSCRIPT or not feature_type:
        for transcript_id in CM().normalize_transcript_id(feature):
            normalized.append((transcript_id, TRANSCRIPT))

    # if the feature is not known by an alias, return the original name
    if not normalized:
        if feature_type:
            normalized = [(feature, feature_type)]
        else:
            raise ValueError(f"Could not determine the type of feature of '{feature}'")

    return normalized


def is_contig(feature: str) -> bool:
    """Return True is the given feature is a contig in the current release."""
    try:
        return any([i[1] == CONTIG for i in normalize_feature(feature)])
    except ValueError:
        return False


def is_exon(feature: str) -> bool:
    """Return True is the given feature is a exon in the current release."""
    try:
        return any([i[1] == EXON for i in normalize_feature(feature)])
    except ValueError:
        return False


def is_gene(feature: str) -> bool:
    """Return True is the given feature is a gene in the current release."""
    try:
        return any([i[1] == GENE for i in normalize_feature(feature)])
    except ValueError:
        return False


def is_protein(feature: str) -> bool:
    """Return True is the given feature is a protein in the current release."""
    try:
        return any([i[1] == PROTEIN for i in normalize_feature(feature)])
    except ValueError:
        return False


def is_transcript(feature: str) -> bool:
    """Return True is the given feature is a transcript in the current release."""
    try:
        return any([i[1] == TRANSCRIPT for i in normalize_feature(feature)])
    except ValueError:
        return False


# --------------------------------------
#  Contig functions
# --------------------------------------
@lru_cache()
def contig_ids(
    feature: str, feature_type: Optional[str] = None, best_only: bool = False
) -> List[str]:
    all_contig_ids = []

    for feature, feature_type in normalize_feature(feature, feature_type):
        if is_ensembl_id(feature):
            all_contig_ids.extend(_get_contig_ids_by_id(feature, feature_type))
        else:
            all_contig_ids.extend(_get_contig_ids_by_name(feature, feature_type, best_only))

    return sorted(filter(None, all_contig_ids))


def _get_contig_ids_by_id(feature: str, feature_type: str) -> List[str]:
    if feature_type == CDS or feature_type == TRANSCRIPT:
        return CM().contig_ids_of_transcript_id(feature)
    elif feature_type == EXON:
        return CM().contig_ids_of_exon_id(feature)
    elif feature_type == GENE:
        return CM().contig_ids_of_gene_id(feature)
    elif feature_type == PROTEIN:
        return CM().contig_ids_of_protein_id(feature)
    else:
        raise ValueError(f"Cannot get contig IDs from (ID={feature}, type={feature_type})")


def _get_contig_ids_by_name(feature: str, feature_type: str, best_only: bool = False) -> List[str]:
    if feature_type == CDS or feature_type == TRANSCRIPT:
        results = []
        for transcript_id in transcript_ids(feature, TRANSCRIPT, best_only):
            results.extend(_get_contig_ids_by_id(transcript_id, TRANSCRIPT))
        return results
    elif feature_type == CONTIG:
        return [feature]
    elif feature_type == GENE:
        results = []
        for gene_id in gene_ids(feature, GENE):
            results.extend(_get_contig_ids_by_id(gene_id, GENE))
        return results
    else:
        raise ValueError(f"Cannot get contig IDs from (name={feature}, type={feature_type})")


# --------------------------------------
#  Exon functions
# --------------------------------------
@lru_cache()
def get_exons(
    feature: str, feature_type: Optional[str] = None, best_only: bool = False
) -> List[pyensembl.Exon]:
    exons = []

    for feature, feature_type in normalize_feature(feature, feature_type):
        for exon_id in exon_ids(feature, feature_type, best_only):
            exons.append(_query(exon_id, CM().ensembl.exon_by_id))

    return exons


@lru_cache()
def exon_ids(
    feature: str, feature_type: Optional[str] = None, best_only: bool = False
) -> List[str]:
    all_exon_ids = []

    for feature, feature_type in normalize_feature(feature, feature_type):
        if is_ensembl_id(feature):
            all_exon_ids.extend(_get_exon_ids_by_id(feature, feature_type, best_only))
        else:
            all_exon_ids.extend(_get_exon_ids_by_name(feature, feature_type))

    return sorted(filter(None, all_exon_ids))


def _get_exon_ids_by_id(feature: str, feature_type: str, best_only: bool = False) -> List[str]:
    if feature_type == CDS or feature_type == TRANSCRIPT:
        return CM().exon_ids_of_transcript_id(feature)
    elif feature_type == EXON:
        return [feature]
    elif feature_type == GENE:
        return CM().exon_ids_of_gene_id(feature)
    elif feature_type == PROTEIN:
        all_exon_ids = []
        for transcript_id in transcript_ids(feature, PROTEIN, best_only):
            all_exon_ids.extend(_get_exon_ids_by_id(transcript_id, TRANSCRIPT, best_only))
        return all_exon_ids
    else:
        raise ValueError(f"Cannot get exon IDs from (ID={feature}, type={feature_type})")


def _get_exon_ids_by_name(feature: str, feature_type: str) -> List[str]:
    if feature_type == CDS or feature_type == TRANSCRIPT:
        return CM().exon_ids_of_transcript_name(feature)
    elif feature_type == CONTIG:
        return CM().exon_ids_of_contig_id(feature)
    elif feature_type == GENE:
        return CM().exon_ids_of_gene_name(feature)
    else:
        raise ValueError(f"Cannot get exon IDs from (name={feature}, type={feature_type})")


# --------------------------------------
#  Gene functions
# --------------------------------------
@lru_cache()
def get_genes(
    feature: str, feature_type: Optional[str] = None, best_only: bool = False
) -> List[pyensembl.Gene]:
    genes = []

    for feature, feature_type in normalize_feature(feature, feature_type):
        for gene_id in gene_ids(feature, feature_type):
            genes.append(_query(gene_id, CM().ensembl.gene_by_id))

    return genes


@lru_cache()
def gene_names(
    feature: str, feature_type: Optional[str] = None, best_only: bool = False
) -> List[str]:
    all_gene_names = []

    for feature, feature_type in normalize_feature(feature, feature_type):
        for gene_id in gene_ids(feature, feature_type):
            all_gene_names.append(CM().ensembl.gene_name_of_gene_id(gene_id))

    return all_gene_names


@lru_cache()
def gene_ids(feature: str, feature_type: Optional[str] = None) -> List[str]:
    all_gene_ids = []

    for feature, feature_type in normalize_feature(feature, feature_type):
        if is_ensembl_id(feature):
            all_gene_ids.extend(_get_gene_ids_by_id(feature, feature_type))
        else:
            all_gene_ids.extend(_get_gene_ids_by_name(feature, feature_type))

    return sorted(filter(None, all_gene_ids))


def _get_gene_ids_by_id(feature: str, feature_type: str) -> List[str]:
    if feature_type == CDS or feature_type == TRANSCRIPT:
        gene_name = _query(feature, CM().ensembl.gene_name_of_transcript_id)
        return _gene_name_to_id(gene_name) if gene_name else []
    elif feature_type == EXON:
        gene_name = _query(feature, CM().ensembl.gene_name_of_exon_id)
        return _gene_name_to_id(gene_name) if gene_name else []
    elif feature_type == GENE:
        return [feature]
    elif feature_type == PROTEIN:
        return [_query(feature, CM().ensembl.gene_id_of_protein_id)]
    else:
        raise ValueError(f"Cannot get gene IDs from (ID={feature}, type={feature_type})")


def _get_gene_ids_by_name(feature: str, feature_type: str) -> List[str]:
    if feature_type == CDS or feature_type == TRANSCRIPT:
        gene_name = _query(feature, CM().ensembl.gene_name_of_transcript_name)
        return _gene_name_to_id(gene_name) if gene_name else []
    elif feature_type == CONTIG:
        return _query(feature, CM().ensembl.all_gene_ids)
    elif feature_type == GENE:
        return _query(feature, CM().ensembl.gene_ids_of_gene_name)
    else:
        raise ValueError(f"Cannot get gene IDs from (name={feature}, type={feature_type})")


def _gene_name_to_id(gene_name: str) -> List[str]:
    return CM().ensembl.gene_ids_of_gene_name(gene_name)


# --------------------------------------
#  Protein functions
# --------------------------------------
@lru_cache()
def protein_ids(
    feature: str, feature_type: Optional[str] = None, best_only: bool = False
) -> List[str]:
    all_protein_ids = []

    for feature, feature_type in normalize_feature(feature, feature_type):
        for transcript in get_transcripts(feature, feature_type, best_only):
            if transcript.protein_id:
                all_protein_ids.append(transcript.protein_id)

    return sorted(filter(None, all_protein_ids))


# --------------------------------------
#  Transcript functions
# --------------------------------------
@lru_cache()
def get_transcripts(
    feature: str, feature_type: Optional[str] = None, best_only: bool = False
) -> List[pyensembl.Transcript]:
    transcripts = []

    for feature, feature_type in normalize_feature(feature, feature_type):
        for transcript_id in transcript_ids(feature, feature_type, best_only):
            transcripts.append(_query(transcript_id, CM().ensembl.transcript_by_id))

    return transcripts


@lru_cache()
def transcript_names(
    feature: str, feature_type: Optional[str] = None, best_only: bool = False
) -> List[str]:
    all_transcript_names = []

    for feature, feature_type in normalize_feature(feature, feature_type):
        for transcript_id in transcript_ids(feature, feature_type, best_only):
            all_transcript_names.append(
                CM().ensembl.transcript_name_of_transcript_id(transcript_id)
            )

    return all_transcript_names


@lru_cache()
def transcript_ids(
    feature: str, feature_type: Optional[str] = None, best_only: bool = False
) -> List[str]:
    all_transcript_ids = []

    for feature, feature_type in normalize_feature(feature, feature_type):
        if is_ensembl_id(feature):
            all_transcript_ids.extend(_get_transcript_ids_by_id(feature, feature_type, best_only))
        else:
            all_transcript_ids.extend(_get_transcript_ids_by_name(feature, feature_type))

    if best_only:
        all_transcript_ids = [i for i in all_transcript_ids if BestTranscript.is_best(i)]

    return sorted(filter(None, all_transcript_ids))


def _get_transcript_ids_with_exon(feature: str, best_only: bool = False) -> List[str]:
    # NOTE: with `pyensembl==1.8.5` calling `transcript_ids_of_exon_ids` does not
    # match anything. As a workaround, we can map the exon to its gene then return
    # all transcripts of that gene that contain the exon.
    all_transcript_ids = []
    exon = _query(feature, CM().ensembl.exon_by_id)
    for transcript in get_transcripts(exon.gene_id, GENE, best_only):
        if feature in [i.exon_id for i in transcript.exons]:
            all_transcript_ids.append(transcript.transcript_id)

    return all_transcript_ids


def _get_transcript_ids_by_id(
    feature: str, feature_type: str, best_only: bool = False
) -> List[str]:
    if feature_type == CDS or feature_type == TRANSCRIPT:
        return [feature]
    elif feature_type == EXON:
        return _get_transcript_ids_with_exon(feature, best_only)
    elif feature_type == GENE:
        return _query(feature, CM().ensembl.transcript_ids_of_gene_id)
    elif feature_type == PROTEIN:
        return [_query(feature, CM().ensembl.transcript_id_of_protein_id)]
    else:
        raise ValueError(f"Cannot get transcript IDs from (ID={feature}, type={feature_type})")


def _get_transcript_ids_by_name(feature: str, feature_type: str) -> List[str]:
    if feature_type == CDS or feature_type == TRANSCRIPT:
        return _query(feature, CM().ensembl.transcript_ids_of_transcript_name)
    elif feature_type == CONTIG:
        return _query(feature, CM().ensembl.all_transcript_ids)
    elif feature_type == GENE:
        return _query(feature, CM().ensembl.transcript_ids_of_gene_name)
    else:
        raise ValueError(f"Cannot get transcript IDs from (name={feature}, type={feature_type})")


# --------------------------------------
#  Query functions
# --------------------------------------
def _query(feature: str, func: Callable) -> Any:
    try:
        return func(feature)
    except ValueError:
        return None
