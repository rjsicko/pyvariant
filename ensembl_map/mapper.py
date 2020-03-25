import logging

from .cache import Ensembl
from .convert import get_map_function
from .features import get_parse_function
from .transcript import get_transcripts
from .util import assert_valid_position


ENSEMBL = Ensembl()


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
    #     ENSEMBL.data.exon_by_id(feature).start,
    #     ENSEMBL.data.exon_by_id(feature).end.end,
    #     "exon",
    #     "transcript",
    # )

    # As a workaround, we can map the exon to its gene then return all transcripts of
    # that gene that contain the exon:

    def contains_exon(transcript_id):
        """Check if the given transcript contains the exon."""
        transcript = ENSEMBL.data.transcript_by_id(transcript_id)
        if feature in [i.exon_id for i in transcript.exons]:
            return True
        else:
            return False

    result = []

    try:
        exon = ENSEMBL.data.exon_by_id(feature)
    except TypeError:
        logging.error(f"No exon '{feature}' found")
        return result

    for pos in gene_to_transcript(exon.gene_id, exon.start, exon.end):
        if contains_exon(pos.transcript_id):
            result.append(pos)

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
    for transcript in transcripts:
        retval = None
        try:
            retval = map_func(transcript, start, end)
            logging.debug(f"Got: {retval}")
        except Exception as exc:
            logging.debug(exc)
            continue
        ret_obj = parse_func(transcript, *retval)
        logging.debug(f"Parsed {transcript, retval} to {ret_obj}")
        result.append(ret_obj)

    return result
