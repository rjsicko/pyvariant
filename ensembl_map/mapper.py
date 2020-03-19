import logging

from .cache import Cache
from .convert import (
    cpos_to_ppos,
    cpos_to_tpos,
    epos_to_gpos,
    gpos_to_tpos,
    ppos_to_cpos,
    gpos_to_epos,
    tpos_to_cpos,
    tpos_to_gpos,
)
from .features import CDS, Exon, Gene, Protein, Transcript
from .util import is_ensembl_id


def cds_to_exon(feature, start, end=None):
    """Map CDS coordinates to exon coordinates."""
    result = []
    for pos in cds_to_transcript(feature, start, end):
        result.extend(transcript_to_exon(pos.transcript_id, pos.start, pos.end))
    return result


def cds_to_gene(feature, start, end=None):
    """Map CDS coordinates to gene coordinates."""
    result = []
    for pos in cds_to_transcript(feature, start, end):
        result.extend(transcript_to_gene(pos.transcript_id, pos.start, pos.end))
    return result


def cds_to_protein(feature, start, end=None):
    """Map CDS coordinates to protein coordinates."""
    return _map(
        feature,
        start,
        end,
        "protein_id",
        None,
        Cache.get_cache().transcript_ids_of_transcript_name,
        cpos_to_ppos,
        Protein.load,
    )


def cds_to_transcript(feature, start, end=None):
    """Map CDS coordinates to transcript coordinates."""
    return _map(
        feature,
        start,
        end,
        "transcript_id",
        None,
        Cache.get_cache().transcript_ids_of_transcript_name,
        cpos_to_tpos,
        Transcript.load,
    )


def exon_to_cds(feature):
    """Map an exon to CDS coordinates."""
    result = []
    for pos in exon_to_transcript(feature):
        result.extend(transcript_to_cds(pos.transcript_id, pos.start, pos.end))
    return result


def exon_to_gene(feature):
    """Map an exon to gene coordinates."""
    result = []
    for pos in exon_to_transcript(feature):
        result.extend(transcript_to_gene(pos.transcript_id, pos.start, pos.end))
    return result


def exon_to_protein(feature):
    """Map an exon to gene coordinates."""
    result = []
    for pos in exon_to_cds(feature):
        result.extend(cds_to_protein(pos.transcript_id, pos.start, pos.end))
    return result


def exon_to_transcript(feature):
    """Map an exon to transcript coordinates."""
    # NOTE: with `pyensembl==1.8.5` calling `transcript_ids_of_exon_ids` does not
    # match anything:
    # return _map(
    #     feature,
    #     Cache.get_cache().exon_by_id(feature).start,
    #     Cache.get_cache().exon_by_id(feature).end,
    #     "transcript_id",
    #     Cache.get_cache().transcript_ids_of_exon_id,
    #     None,
    #     gpos_to_tpos,
    # )

    # As a workaround, we can map the exon to its gene then return all transcripts of
    # that gene that contain the exon:
    transcripts_with_exon = []
    exon = Cache.get_cache().exon_by_id(feature)
    transcript_ids = Cache.get_cache().transcript_ids_of_gene_id(exon.gene_id)
    for tid in transcript_ids:
        transcript = Cache.get_cache().transcript_by_id(tid)
        if feature in [i.exon_id for i in transcript.exons]:
            transcripts_with_exon.append(tid)

    result = []
    for tid in transcripts_with_exon:
        result.extend(
            _map(
                tid,
                exon.start,
                exon.end,
                "transcript_id",
                None,
                None,
                gpos_to_tpos,
                Transcript.load,
            )
        )
    return result


def gene_to_cds(feature, start, end=None):
    """Map gene coordinates to CDS coordinates."""
    result = []
    for pos in gene_to_transcript(feature, start, end):
        result.extend(transcript_to_cds(pos.transcript_id, pos.start, pos.end))
    return result


def gene_to_exon(feature, start, end=None):
    """Map gene coordinates to exon coordinates."""
    return _map(
        feature,
        start,
        end,
        "gene_id",
        Cache.get_cache().transcript_ids_of_gene_id,
        Cache.get_cache().transcript_ids_of_gene_name,
        gpos_to_epos,
        Exon.load,
    )


def gene_to_protein(feature, start, end=None):
    """Map gene coordinates to protein coordinates."""
    result = []
    for pos in gene_to_cds(feature, start, end):
        result.extend(cds_to_protein(pos.transcript_id, pos.start, pos.end))
    return result


def gene_to_transcript(feature, start, end=None):
    """Map gene coordinates to transcript coordinates."""
    return _map(
        feature,
        start,
        end,
        "transcript_id",
        Cache.get_cache().transcript_ids_of_gene_id,
        Cache.get_cache().transcript_ids_of_gene_name,
        gpos_to_tpos,
        Transcript.load,
    )


def protein_to_cds(feature, start, end=None):
    """Map protein coordinates to CDS coordinates."""
    cds = _map(
        feature,
        start,
        end,
        "transcript_id",
        Cache.get_cache().transcript_id_of_protein_id,
        None,
        ppos_to_cpos,
        CDS.load,
    )
    # add 2 to the end position to get the end of the codon
    for i in cds:
        i.end += 2
    return cds


def protein_to_exon(feature, start, end=None):
    """Map protein coordinates to exon coordinates."""
    result = []
    for pos in protein_to_transcript(feature, start, end):
        result.extend(transcript_to_exon(pos.transcript_id, pos.start, pos.end))
    return result


def protein_to_gene(feature, start, end=None):
    """Map protein coordinates to gene coordinates."""
    result = []
    for pos in protein_to_transcript(feature, start, end):
        result.extend(transcript_to_gene(pos.transcript_id, pos.start, pos.end))
    return result


def protein_to_transcript(feature, start, end=None):
    """Map protein coordinates to transcript coordinates."""
    result = []
    for pos in protein_to_cds(feature, start, end):
        result.extend(cds_to_transcript(pos.transcript_id, pos.start, pos.end))
    return result


def transcript_to_cds(feature, start, end=None):
    """Map transcript coordinates to CDS coordinates."""
    return _map(
        feature,
        start,
        end,
        "transcript_id",
        None,
        Cache.get_cache().transcript_ids_of_transcript_name,
        tpos_to_cpos,
        CDS.load,
    )


def transcript_to_exon(feature, start, end=None):
    """Map transcript coordinates to exon coordinates."""
    result = []
    for pos in transcript_to_gene(feature, start, end):
        result.extend(gene_to_exon(pos.gene_id, pos.start, pos.end))
    return result


def transcript_to_gene(feature, start, end=None):
    """Map transcript coordinates to gene coordinates."""
    return _map(
        feature,
        start,
        end,
        "gene_id",
        None,
        Cache.get_cache().transcript_ids_of_transcript_name,
        tpos_to_gpos,
        Gene.load,
    )


def transcript_to_protein(feature, start, end=None):
    """Map transcript coordinates to protein coordinates."""
    result = []
    for pos in transcript_to_cds(feature, start, end):
        result.extend(cds_to_protein(pos.transcript_id, pos.start, pos.end))
    return result


def _map(feature, start, end, return_id, tids_by_id, tids_by_name, convert, load):
    """
    Template function for mapping a feature to the associated transcript(s), then 
    converting the given coordinates.

    Args:
        feature (str): feature name or Ensembl ID
        start (int): first position relative to `feature`
        end (int or None): second position relative to `feature`
        return_id (str): `Transcript` attribute to return
        tids_by_id: function converts an Ensembl ID to a list of Ensembl transcript IDs
        tids_by_name: function converts a name to a list of Ensembl transcript IDs
        convert: function maps the given coordinate type to the desired type
        load: parse the returned values into a python object
        
    Returns:
        list: of converted coordinate (key, start, end, strand)
    """
    transcript_ids = []
    result = []

    logging.debug(
        "From: "
        f"{feature}, "
        f"{start}, "
        f"{end}, "
        f"{return_id}, "
        f"{tids_by_id}, "
        f"{tids_by_name}, "
        f"{convert}"
    )

    if start < 1:
        raise ValueError(f"Start must be >= 1 ({start})")
    if end is not None and end < 1:
        raise ValueError(f"End must be >= 1 ({end})")

    # map feature to Ensembl transcript ID
    if is_ensembl_id(feature):
        if tids_by_id:
            logging.debug(f"[{tids_by_id}, '{feature}']: {transcript_ids}")
            transcript_ids = tids_by_id(feature)
        else:
            transcript_ids = feature
    else:
        if tids_by_name:
            logging.debug(f"[{tids_by_name}, '{feature}']: {transcript_ids}")
            transcript_ids = tids_by_name(feature)

    if not transcript_ids:
        raise ValueError(f"'{feature}' did not match any transcript ID")

    # coerce to a list
    if isinstance(transcript_ids, str):
        transcript_ids = [transcript_ids]

    # map coordinates by transcript ID
    for tid in sorted(transcript_ids):
        args = []
        ret_val1 = None
        ret_val2 = None
        tobj = Cache.get_cache().transcript_by_id(tid)
        try:
            if start is not None:
                ret_val1 = convert(tobj, start)
            if end is not None:
                ret_val2 = convert(tobj, end)
            logging.debug(f"Got: {ret_val1} {ret_val2}")
        except Exception as exc:
            logging.debug(exc)
            continue

        if isinstance(ret_val1, int):
            args = [ret_val1, ret_val2]
            if ret_val2 is not None:
                args = [ret_val1, ret_val2]
            else:
                args = [ret_val1, ret_val1]
        else:
            args = ret_val1

        ret_obj = load(tobj, *args)
        logging.debug(f"Load: {ret_obj}")
        result.append(ret_obj)

    return result
