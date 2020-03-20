from .cache import Cache
from .util import is_ensembl_id


def get_transcripts(feature, feature_type):
    transcripts = []
    for transcript_id in _get_transcript_ids(feature, feature_type):
        transcripts.append(Cache.get_cache().transcript_by_id(transcript_id))
    return transcripts


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
