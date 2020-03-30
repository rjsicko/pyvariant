from .cache import Ensembl
from .util import is_ensembl_id


ENSEMBL = Ensembl()


def get_transcripts(feature, feature_type):
    transcripts = []
    for transcript_id in _get_transcript_ids(feature, feature_type):
        transcripts.append(_query(transcript_id, "transcript", ENSEMBL.data.transcript_by_id))
    return transcripts


def _get_transcript_ids(feature, feature_type):
    if is_ensembl_id(feature):
        transcript_ids = _get_transcript_ids_by_id(feature, feature_type)
    else:
        transcript_ids = _get_transcript_ids_by_name(feature, feature_type)

    if not isinstance(transcript_ids, list):
        transcript_ids = [transcript_ids]

    return sorted(transcript_ids)


def _get_transcript_ids_with_exon(feature):
    # NOTE: with `pyensembl==1.8.5` calling `transcript_ids_of_exon_ids` does not
    # match anything. As a workaround, we can map the exon to its gene then return
    # all transcripts of that gene that contain the exon.
    transcript_ids = []
    exon = _query(feature, "exon", ENSEMBL.data.exon_by_id)
    for transcript in get_transcripts(exon.gene_id, "gene"):
        if feature in [i.exon_id for i in transcript.exons]:
            transcript_ids.append(transcript.transcript_id)
    return transcript_ids


def _get_transcript_ids_by_id(feature, feature_type):
    if feature_type == "cds" or feature_type == "transcript":
        return feature
    elif feature_type == "exon":
        return _get_transcript_ids_with_exon(feature)
    elif feature_type == "gene":
        return _query(feature, feature_type, ENSEMBL.data.transcript_ids_of_gene_id)
    elif feature_type == "protein":
        return _query(feature, feature_type, ENSEMBL.data.transcript_id_of_protein_id)
    else:
        raise TypeError(f"Cannot get transcript IDs from {feature_type}")


def _get_transcript_ids_by_name(feature, feature_type):
    if feature_type == "cds" or feature_type == "transcript":
        return _query(feature, "transcript", ENSEMBL.data.transcript_ids_of_transcript_name)
    elif feature_type == "gene":
        return _query(feature, feature_type, ENSEMBL.data.transcript_ids_of_gene_name)
    else:
        raise TypeError(f"Cannot get transcript IDs from {feature_type}")


def _query(feature, feature_type, func):
    try:
        return func(feature)
    except ValueError:
        raise ValueError(f"No match for {feature_type} '{feature}'")
