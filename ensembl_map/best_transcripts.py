import csv
from typing import Dict, Optional

from logzero import logger


class BestTranscript:
    """Helper class for gene-to-best-transcript mappings."""

    _data: Dict[str, str] = {}
    _loaded = False

    @classmethod
    def load(cls, path: str = ""):
        """Load one-to-one mappings of Ensembl gene IDs to transcript IDs.

        Expects a tab-seperated file like:
            Ensembl_Gene_ID Ensembl_Transcript_ID
            ENSG00000276626 ENST00000612820
            ENSG00000201317 ENST00000364447
            ENSG00000200823 ENST00000363953
            ...

        Args:
            path (str): path to the file

        Example:
            >>> BestTranscript.load()
            >>> BestTranscript.get("ENSG00000207588")
            'ENST00000384856'
        """
        cls._data = {}
        cls._loaded = False

        if path:
            logger.debug(f"Loading best transcripts from '{path}'")
            with open(path, "r") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    cls._data[row["Ensembl_Gene_ID"]] = row["Ensembl_Transcript_ID"]
        else:
            logger.warning("No best transcript file")

        cls._loaded = True

    @classmethod
    def get(cls, gene_id: str) -> Optional[str]:
        """Return the best transcript for the given gene."""
        if not cls._loaded:
            cls.load()

        return cls._data.get(gene_id, None)

    @classmethod
    def is_best(cls, transcript_id: str) -> bool:
        """The given transcript is in the list of best transcripts."""
        if not cls._loaded:
            cls.load()

        return transcript_id in cls._data.values()


def get_best_transcript(gene: str) -> Optional[str]:
    """Return the best transcript for the given gene."""
    return BestTranscript.get(gene)


def is_best_transcript(transcript: str) -> bool:
    """Return True if the given transcript is in the list of best transcripts."""
    return BestTranscript.is_best(transcript)
