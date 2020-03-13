import logging

from .cache import data
from .convert import (
    cpos_to_ppos,
    cpos_to_tpos,
    epos_to_tpos,
    gpos_to_tpos,
    ppos_to_cpos,
    tpos_to_epos,
    tpos_to_cpos,
    tpos_to_gpos,
)
from .util import is_ensembl_id


def cds_to_exon(feature, position, end=None):
    """Map CDS coordinates to exon position."""
    result = []
    for pos in cds_to_transcript(feature, position):
        result.extend(transcript_to_exon(*pos[:3]))

    return result


def cds_to_gene(feature, position, end=None):
    """Map CDS coordinates to gene coordinates."""
    result = []
    for pos in cds_to_transcript(feature, position):
        result.extend(transcript_to_gene(*pos[:3]))

    return result


def cds_to_protein(feature, position, end=None):
    """Map CDS coordinates to protein coordinates."""
    return _map(
        feature,
        position,
        None,
        "protein_id",
        None,
        data.transcript_ids_of_transcript_name,
        cpos_to_ppos,
    )


def cds_to_transcript(feature, position, end=None):
    """Map CDS coordinates to transcript coordinates."""
    return _map(
        feature,
        position,
        None,
        "transcript_id",
        None,
        data.transcript_ids_of_transcript_name,
        cpos_to_tpos,
    )


def exon_to_cds(feature, position):
    """Map exon position to CDS coordinates."""
    result = []
    for pos in exon_to_transcript(feature, position):
        result.extend(transcript_to_cds(*pos[:3]))

    return result


def exon_to_gene(feature, position):
    """Map exon position to gene coordinates."""
    result = []
    for pos in exon_to_transcript(feature, position):
        result.extend(transcript_to_gene(*pos[:3]))

    return result


def exon_to_protein(feature, position):
    """Map exon position to gene coordinates."""
    result = []
    for pos in exon_to_cds(feature, position):
        result.extend(cds_to_protein(*pos[:3]))

    return result


def exon_to_transcript(feature, position):
    """Map exon position to transcript coordinates."""
    return _map(
        feature,
        position,
        None,
        "transcript_id",
        data.transcript_ids_of_exon_id,
        None,
        epos_to_tpos,
    )


def gene_to_cds(feature, position, end=None):
    """Map gene coordinates to CDS coordinates."""
    result = []
    for pos in gene_to_transcript(feature, position):
        result.extend(transcript_to_cds(*pos[:3]))

    return result


def gene_to_exon(feature, position, end=None):
    """Map gene coordinates to exon position."""
    raise NotImplementedError


def gene_to_protein(feature, position, end=None):
    """Map gene coordinates to protein coordinates."""
    result = []
    for pos in gene_to_cds(feature, position, end=None):
        result.extend(cds_to_protein(*pos[:3]))

    return result


def gene_to_transcript(feature, position, end=None):
    """Map gene coordinates to transcript coordinates."""
    return _map(
        feature,
        position,
        end,
        "transcript_id",
        data.transcript_ids_of_gene_id,
        data.transcript_ids_of_gene_name,
        gpos_to_tpos,
    )


def protein_to_cds(feature, position, end=None):
    """Map protein coordinates to CDS coordinates."""
    return _map(
        feature,
        position,
        end,
        "transcript_id",
        data.transcript_id_of_protein_id,
        None,
        ppos_to_cpos,
    )


def protein_to_exon(feature, position, end=None):
    """Map protein coordinates to exon position."""
    result = []
    for pos in protein_to_transcript(feature, position):
        result.extend(transcript_to_exon(*pos[:3]))

    return result


def protein_to_gene(feature, position, end=None):
    """Map protein coordinates to gene coordinates."""
    result = []
    for pos in protein_to_transcript(feature, position, end=None):
        result.extend(transcript_to_gene(*pos[:3]))

    return result


def protein_to_transcript(feature, position, end=None):
    """Map protein coordinates to transcript coordinates."""
    result = []
    for pos in protein_to_cds(feature, position, end=None):
        result.extend(cds_to_gene(*pos[:3]))

    return result


def transcript_to_cds(feature, position, end=None):
    """Map transcript coordinates to CDS coordinates."""
    return _map(
        feature,
        position,
        end,
        "transcript_id",
        None,
        data.transcript_ids_of_transcript_name,
        tpos_to_cpos,
    )


def transcript_to_exon(feature, position, end=None):
    """Map transcript coordinates to exon position."""
    return _map(
        feature,
        position,
        end,
        "gene_id",
        None,
        data.transcript_ids_of_transcript_name,
        tpos_to_epos,
    )


def transcript_to_gene(feature, position, end=None):
    """Map transcript coordinates to gene coordinates."""
    return _map(
        feature,
        position,
        end,
        "gene_id",
        None,
        data.transcript_ids_of_transcript_name,
        tpos_to_gpos,
    )


def transcript_to_protein(feature, position, end=None):
    """Map transcript coordinates to protein coordinates."""
    result = []
    for pos in transcript_to_cds(feature, position, end=None):
        result.extend(cds_to_protein(*pos[:3]))

    return result


def _map(feature, position, end, return_id, tids_by_id, tids_by_name, convert):
    """
    Template function for mapping a feature to the associated transcript(s), then 
    converting the given coordinates.

    Args:
        feature (str): feature name or Ensembl ID
        position (int): first position relative to `feature`
        end (int or None): second position relative to `feature`
        return_id (str): `Transcript` attribute to return
        tids_by_id: function converts an Ensembl ID to a list of Ensembl transcript IDs
        tids_by_name: function converts a name to a list of Ensembl transcript IDs
        convert: function maps the given coordinate type to the desired type
        
    Returns:
        list: of converted coordinate (key, position, end, strand)
    """
    transcript_ids = []
    result = []

    assert feature is not None, f"Feature must not be 'None'"
    assert position >= 1, f"Position is not >= 1 ({position})"
    if end:
        assert end >= 1, f"End is not >= 1 ({position})"

    # map feature to Ensembl transcript ID
    if is_ensembl_id(feature):
        if tids_by_id:
            transcript_ids = tids_by_id(feature)
            logging.debug(f"From '{feature}' got transcript IDs: {transcript_ids}")
        else:
            transcript_ids = [feature]
    else:
        if tids_by_name:
            transcript_ids = tids_by_name(feature)
            logging.debug(f"From '{feature}' got transcript IDs: {transcript_ids}")

    if not transcript_ids:
        raise ValueError(f"Query did not return any transcript IDs")

    # map coordinates by transcript ID
    for tid in sorted(transcript_ids):
        tobj = data.transcript_by_id(tid)
        try:
            pos1 = convert(tobj, position)
            if end:
                pos2 = convert(tobj, end)
            else:
                pos2 = pos1
        except Exception as exc:
            logging.debug(exc)
            continue
        # return coordinates relative to the positive strand
        if pos2 < pos1:
            pos1, pos2 = pos2, pos1
        result.append((getattr(tobj, return_id), pos1, pos2, tobj.strand))

    return result
