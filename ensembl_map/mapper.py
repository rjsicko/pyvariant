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
    return _map(feature, start, end, "cds", "protein")


def cds_to_transcript(feature, start, end=None):
    """Map CDS coordinates to transcript coordinates."""
    return _map(feature, start, end, "cds", "transcript")


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
    #     Cache.get_cache().exon_by_id(feature).end.end,
    #     "exon",
    #     "transcript",
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
        result.extend(_map(tid, exon.start, exon.end, "exon", "transcript"))
    return result


def gene_to_cds(feature, start, end=None):
    """Map gene coordinates to CDS coordinates."""
    result = []
    for pos in gene_to_transcript(feature, start, end):
        result.extend(transcript_to_cds(pos.transcript_id, pos.start, pos.end))
    return result


def gene_to_exon(feature, start, end=None):
    """Map gene coordinates to exon coordinates."""
    return _map(feature, start, end, "gene", "exon")


def gene_to_protein(feature, start, end=None):
    """Map gene coordinates to protein coordinates."""
    result = []
    for pos in gene_to_cds(feature, start, end):
        result.extend(cds_to_protein(pos.transcript_id, pos.start, pos.end))
    return result


def gene_to_transcript(feature, start, end=None):
    """Map gene coordinates to transcript coordinates."""
    return _map(feature, start, end, "gene", "transcript")


def protein_to_cds(feature, start, end=None):
    """Map protein coordinates to CDS coordinates."""
    cds = _map(feature, start, end, "protein", "cds")
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
    return _map(feature, start, end, "transcript", "cds")


def transcript_to_exon(feature, start, end=None):
    """Map transcript coordinates to exon coordinates."""
    result = []
    for pos in transcript_to_gene(feature, start, end):
        result.extend(gene_to_exon(pos.gene_id, pos.start, pos.end))
    return result


def transcript_to_gene(feature, start, end=None):
    """Map transcript coordinates to gene coordinates."""
    return _map(feature, start, end, "transcript", "gene")


def transcript_to_protein(feature, start, end=None):
    """Map transcript coordinates to protein coordinates."""
    result = []
    for pos in transcript_to_cds(feature, start, end):
        result.extend(cds_to_protein(pos.transcript_id, pos.start, pos.end))
    return result


def _map(feature, start, end, from_type, to_type):
    """
    Template function for mapping a feature to the associated transcript(s), then 
    converting the given coordinates.

    Args:
        feature (str): feature name or Ensembl ID
        start (int): first position relative to `feature`
        end (int or None): second position relative to `feature`
        from_type (str): coordinates relative to this type of feature (e.g. 'gene')
        to_type (str): map coordinates to this type of feature (e.g 'transcript')
        
    Returns:
        list: of converted coordinate
    """
    transcript_ids = []
    result = []

    logging.debug(f"Map ({feature}, {start}, {end}) from {from_type} to {to_type}")

    _assert_valid_position(start, end)
    map_func = _get_map_function(from_type, to_type)
    parse_func = _get_parse_function(to_type)
    transcript_ids = _get_transcript_ids(feature, from_type)

    # map coordinates by transcript ID
    for tid in transcript_ids:
        args = []
        ret_val1 = None
        ret_val2 = None
        tobj = Cache.get_cache().transcript_by_id(tid)
        try:
            if start is not None:
                ret_val1 = map_func(tobj, start)
            if end is not None:
                ret_val2 = map_func(tobj, end)
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

        ret_obj = parse_func(tobj, *args)
        logging.debug(f"Parsed: {ret_obj}")
        result.append(ret_obj)

    return result


def _assert_valid_position(start=None, end=None):
    if start is not None and start < 1:
        raise ValueError(f"Start must be >= 1 ({start})")
    if end is not None and end < 1:
        raise ValueError(f"End must be >= 1 ({end})")


def _get_transcript_ids(feature, feature_type):
    if is_ensembl_id(feature):
        transcript_ids = _get_transcript_ids_by_id(feature, feature_type)
    else:
        transcript_ids = _get_transcript_ids_by_name(feature, feature_type)

    if not transcript_ids:
        raise ValueError(f"'{feature}' did not match any transcript ID")

    if not isinstance(transcript_ids, list):
        transcript_ids = [transcript_ids]

    return sorted(transcript_ids)


def _get_transcript_ids_by_id(feature, feature_type):
    if feature_type == "cds":
        return feature
    elif feature_type == "exon":
        x = Cache.get_cache().exon_by_id(feature).gene_id
        return Cache.get_cache().transcript_ids_of_gene_id(x)
    elif feature_type == "gene":
        return Cache.get_cache().transcript_ids_of_gene_id(feature)
    elif feature_type == "protein":
        return Cache.get_cache().transcript_id_of_protein_id(feature)
    elif feature_type == "transcript":
        return feature
    else:
        raise TypeError(f"Could not get transcript IDs for {feature_type}")


def _get_transcript_ids_by_name(feature, feature_type):
    if feature_type == "cds":
        return Cache.get_cache().transcript_ids_of_transcript_name(feature)
    elif feature_type == "exon":
        raise NotImplementedError(f"Cannot get transcript IDs for {feature_type} name")
    elif feature_type == "gene":
        return Cache.get_cache().transcript_ids_of_gene_name(feature)
    elif feature_type == "protein":
        raise NotImplementedError(f"Cannot get transcript IDs for {feature_type} name")
    elif feature_type == "transcript":
        return Cache.get_cache().transcript_ids_of_transcript_name(feature)
    else:
        raise TypeError(f"Could not get transcript IDs for {feature_type}")


def _get_map_function(from_type, to_type):
    if from_type == to_type:
        raise ValueError("`from_type` and `to_type` must be different!")

    elif from_type == "cds" and to_type == "exon":
        raise NotImplementedError(f"Cannot convert {from_type} directly to {to_type}")
    elif from_type == "cds" and to_type == "gene":
        raise NotImplementedError(f"Cannot convert {from_type} directly to {to_type}")
    elif from_type == "cds" and to_type == "protein":
        return cpos_to_ppos
    elif from_type == "cds" and to_type == "transcript":
        return cpos_to_tpos

    elif from_type == "exon" and to_type == "cds":
        raise NotImplementedError(f"Cannot convert {from_type} directly to {to_type}")
    elif from_type == "exon" and to_type == "gene":
        return epos_to_gpos
    elif from_type == "exon" and to_type == "protein":
        raise NotImplementedError(f"Cannot convert {from_type} directly to {to_type}")
    elif from_type == "exon" and to_type == "transcript":
        raise NotImplementedError(f"Cannot convert {from_type} directly to {to_type}")

    elif from_type == "gene" and to_type == "cds":
        raise NotImplementedError(f"Cannot convert {from_type} directly to {to_type}")
    elif from_type == "gene" and to_type == "exon":
        return gpos_to_epos
    elif from_type == "gene" and to_type == "protein":
        raise NotImplementedError(f"Cannot convert {from_type} directly to {to_type}")
    elif from_type == "gene" and to_type == "transcript":
        return gpos_to_tpos

    elif from_type == "protein" and to_type == "cds":
        return ppos_to_cpos
    elif from_type == "protein" and to_type == "exon":
        raise NotImplementedError(f"Cannot convert {from_type} directly to {to_type}")
    elif from_type == "protein" and to_type == "gene":
        raise NotImplementedError(f"Cannot convert {from_type} directly to {to_type}")
    elif from_type == "protein" and to_type == "transcript":
        raise NotImplementedError(f"Cannot convert {from_type} directly to {to_type}")

    elif from_type == "transcript" and to_type == "cds":
        return tpos_to_cpos
    elif from_type == "transcript" and to_type == "exon":
        raise NotImplementedError(f"Cannot convert {from_type} directly to {to_type}")
    elif from_type == "transcript" and to_type == "gene":
        return tpos_to_gpos
    elif from_type == "transcript" and to_type == "protein":
        raise NotImplementedError(f"Cannot convert {from_type} directly to {to_type}")

    else:
        raise ValueError(f"Unrecognized combination, from '{from_type}' to '{to_type}'")


def _get_parse_function(to_type):
    if to_type == "cds":
        return CDS.load
    elif to_type == "exon":
        return Exon.load
    elif to_type == "gene":
        return Gene.load
    elif to_type == "protein":
        return Protein.load
    elif to_type == "transcript":
        return Transcript.load
    else:
        raise TypeError(f"Could not get parse function for {to_type}")
