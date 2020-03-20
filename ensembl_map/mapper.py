import logging

from .cache import Cache
from .convert import get_map_function
from .features import get_parse_function
from .transcript import get_transcripts
from .util import assert_valid_position


def cds_to_exon(feature, start, end=None):
    """Map CDS coordinates to exon coordinates."""
    result = []
    for pos in cds_to_transcript(feature, start, end):
        result.extend(transcript_to_exon(*pos.to_tuple()))
    return result


def cds_to_gene(feature, start, end=None):
    """Map CDS coordinates to gene coordinates."""
    result = []
    for pos in cds_to_transcript(feature, start, end):
        result.extend(transcript_to_gene(*pos.to_tuple()))
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
        result.extend(transcript_to_cds(*pos.to_tuple()))
    return result


def exon_to_gene(feature):
    """Map an exon to gene coordinates."""
    result = []
    for pos in exon_to_transcript(feature):
        result.extend(transcript_to_gene(*pos.to_tuple()))
    return result


def exon_to_protein(feature):
    """Map an exon to gene coordinates."""
    result = []
    for pos in exon_to_cds(feature):
        result.extend(cds_to_protein(*pos.to_tuple()))
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
        result.extend(transcript_to_cds(*pos.to_tuple()))
    return result


def gene_to_exon(feature, start, end=None):
    """Map gene coordinates to exon coordinates."""
    return _map(feature, start, end, "gene", "exon")


def gene_to_protein(feature, start, end=None):
    """Map gene coordinates to protein coordinates."""
    result = []
    for pos in gene_to_cds(feature, start, end):
        result.extend(cds_to_protein(*pos.to_tuple()))
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
        result.extend(transcript_to_exon(*pos.to_tuple()))
    return result


def protein_to_gene(feature, start, end=None):
    """Map protein coordinates to gene coordinates."""
    result = []
    for pos in protein_to_transcript(feature, start, end):
        result.extend(transcript_to_gene(*pos.to_tuple()))
    return result


def protein_to_transcript(feature, start, end=None):
    """Map protein coordinates to transcript coordinates."""
    result = []
    for pos in protein_to_cds(feature, start, end):
        result.extend(cds_to_transcript(*pos.to_tuple()))
    return result


def transcript_to_cds(feature, start, end=None):
    """Map transcript coordinates to CDS coordinates."""
    return _map(feature, start, end, "transcript", "cds")


def transcript_to_exon(feature, start, end=None):
    """Map transcript coordinates to exon coordinates."""
    result = []
    for pos in transcript_to_gene(feature, start, end):
        result.extend(gene_to_exon(*pos.to_tuple()))
    return result


def transcript_to_gene(feature, start, end=None):
    """Map transcript coordinates to gene coordinates."""
    return _map(feature, start, end, "transcript", "gene")


def transcript_to_protein(feature, start, end=None):
    """Map transcript coordinates to protein coordinates."""
    result = []
    for pos in transcript_to_cds(feature, start, end):
        result.extend(cds_to_protein(*pos.to_tuple()))
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
    result = []

    logging.debug(f"Map {from_type} ({feature}, {start}, {end}) to {to_type}")
    assert_valid_position(start, end)
    map_func = get_map_function(from_type, to_type)
    parse_func = get_parse_function(to_type)
    transcripts = get_transcripts(feature, from_type)

    # map coordinates by transcript ID
    for tobj in transcripts:
        args = [tobj]
        ret_val1 = None
        ret_val2 = None
        try:
            if start is not None:
                ret_val1 = map_func(tobj, start)
                logging.debug(f"Got: {ret_val1}")
            if end is not None:
                ret_val2 = map_func(tobj, end)
                logging.debug(f"Got: {ret_val2}")
        except Exception as exc:
            logging.debug(exc)
            continue

        if isinstance(ret_val1, int):
            if ret_val2 is not None:
                args.extend([ret_val1, ret_val2])
            else:
                args.extend([ret_val1, ret_val1])
        else:
            args.extend(ret_val1)

        ret_obj = parse_func(*args)
        logging.debug(f"Parsed {args} to {ret_obj}")
        result.append(ret_obj)

    return result
