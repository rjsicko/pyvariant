import os.path
from ftplib import FTP
from pathlib import Path

import pandas
from gtfparse import read_gtf
from logzero import logger
from pyfaidx import Fasta

from .utils import bgzip, is_bgzipped, strip_version

# default cache directory
DEFAULT_CACHE_DIR = "."

# Ensembl FTP URL
ENSEMBL_FTP_SERVER = "ftp.ensembl.org"

# FASTA directory example: pub/release-100/fasta/homo_sapiens
FASTA_SUBDIR_TEMPLATE = "pub/release-{release}/fasta/{species}/{type}"

# FASTA example: Homo_sapiens.GRCh38.100.dna.toplevel.fa.gz
# the file name format changed after Ensembl release 75
FASTA_FILENAME_TEMPLATE_OLD = {
    "cdna": "{species}.{reference}.{release}.{type}.all.fa.gz",
    "dna": "{species}.{reference}.{release}.{type}.toplevel.fa.gz",
    "ncrna": "{species}.{reference}.{release}.{type}.fa.gz",
    "pep": "{species}.{reference}.{release}.{type}.all.fa.gz",
}
FASTA_FILENAME_TEMPLATE_NEW = {
    "cdna": "{species}.{reference}.{type}.all.fa.gz",
    "dna": "{species}.{reference}.{type}.toplevel.fa.gz",
    "ncrna": "{species}.{reference}.{type}.fa.gz",
    "pep": "{species}.{reference}.{type}.all.fa.gz",
}

# GTF annotation directory example: pub/release-100/gtf/homo_sapiens
GTF_SUBDIR_TEMPLATE = "pub/release-{release}/gtf/{species}"

# GTF annotation file example: Homo_sapiens.GRCh38.100.gtf.gz
GTF_FILENAME_TEMPLATE = "{species}.{reference}.{release}.gtf.gz"


# load GTF into a pandas dataframe
GTF_COLUMN_RENAME = {"seqname": "contig_id"}

GTF_KEEP_COLUMNS = [
    "contig_id",
    "feature",
    "start",
    "end",
    "strand",
    "gene_id",
    "gene_name",
    "transcript_id",
    "transcript_name",
    "exon_id",
    "exon_number",
    "protein_id",
]

GTF_KEEP_FEATURES = ["CDS", "exon", "stop_codon", "transcript"]


class Cache:
    """Class for managing the data files required by this package."""

    def __init__(
        self, species: str, reference: str, release: int, cache_dir: str = DEFAULT_CACHE_DIR
    ):
        self.species = species
        self.reference = reference
        self.release = release
        self.cache_dir = os.path.abspath(cache_dir)

    @property
    def release_cache_dir(self) -> str:
        """Local path to the directory containing the release data files."""
        return os.path.join(self.cache_dir, self.species, self.reference, str(self.release))

    def make_release_cache_dir(self):
        """Create the release cache directory, if it doesn't already exist."""
        Path(self.release_cache_dir).mkdir(exist_ok=True, parents=True)

    def download_all(self):
        """Download required data and cache it.

        NOTE: pyfaidx only supports compressed FASTA in BGZF format. Ensembl FASTA comes in GZ
        format, so we need to compress the data with bgzip.
        """
        self.make_release_cache_dir()

        for fasta, index, downloadf, indexf in [
            (
                self.local_cdna_fasta_filepath,
                self.local_cdna_index_filepath,
                self.download_cdna_fasta,
                self.index_cdna_fasta,
            ),
            (
                self.local_dna_fasta_filepath,
                self.local_dna_index_filepath,
                self.download_dna_fasta,
                self.index_dna_fasta,
            ),
            (
                self.local_pep_fasta_filepath,
                self.local_pep_index_filepath,
                self.download_pep_fasta,
                self.index_pep_fasta,
            ),
            (
                self.local_ncrna_fasta_filepath,
                self.local_ncrna_index_filepath,
                self.download_ncrna_fasta,
                self.index_ncrna_fasta,
            ),
        ]:
            if not os.path.exists(fasta):
                downloadf()
            if not is_bgzipped(fasta):
                bgzip(fasta)
            if not os.path.exists(index):
                indexf()

        if not os.path.exists(self.local_gtf_cache_filepath):
            if not os.path.exists(self.local_gtf_filepath):
                self.download_gtf()
            df = self.normalize_gtf()
            self.cache_gtf(df)

    # ---------------------------------------------------------------------------------------------
    # FASTA files
    # ---------------------------------------------------------------------------------------------
    @property
    def remote_cdna_fasta_subdir(self) -> str:
        """Path to the directory on the FTP server that contains the cDNA FASTA file."""
        return self._remote_fasta_subdir("cdna")

    @property
    def remote_dna_fasta_subdir(self) -> str:
        """Path to the directory on the FTP server that contains the DNA FASTA file."""
        return self._remote_fasta_subdir("dna")

    @property
    def remote_ncrna_fasta_subdir(self) -> str:
        """Path to the directory on the FTP server that contains the ncRNA FASTA file."""
        return self._remote_fasta_subdir("ncrna")

    @property
    def remote_pep_fasta_subdir(self) -> str:
        """Path to the directory on the FTP server that contains the peptide FASTA file."""
        return self._remote_fasta_subdir("pep")

    def _remote_fasta_subdir(self, fasta_type: str) -> str:
        return FASTA_SUBDIR_TEMPLATE.format(
            release=self.release, species=self.species, type=fasta_type
        )

    @property
    def remote_cdna_fasta_filename(self) -> str:
        """Name of the cDNA FASTA file on the FTP server."""
        return self._remote_fasta_filename("cdna")

    @property
    def remote_dna_fasta_filename(self) -> str:
        """Name of the DNA FASTA file on the FTP server."""
        return self._remote_fasta_filename("dna")

    @property
    def remote_ncrna_fasta_filename(self) -> str:
        """Name of the ncRNA FASTA file on the FTP server."""
        return self._remote_fasta_filename("ncrna")

    @property
    def remote_pep_fasta_filename(self) -> str:
        """Name of the peptide FASTA file on the FTP server."""
        return self._remote_fasta_filename("pep")

    def _remote_fasta_filename(self, fasta_type: str) -> str:
        if self.release <= 75:
            template = FASTA_FILENAME_TEMPLATE_OLD[fasta_type]
            return template.format(
                species=self.species.capitalize(),
                reference=self.reference,
                release=self.release,
                type=fasta_type,
            )
        else:
            template = FASTA_FILENAME_TEMPLATE_NEW[fasta_type]
            return template.format(
                species=self.species.capitalize(), reference=self.reference, type=fasta_type
            )

    @property
    def local_cdna_fasta_filepath(self) -> str:
        """Local name of the cDNA FASTA file."""
        return self._local_fasta_filepath(self.remote_cdna_fasta_filename)

    @property
    def local_dna_fasta_filepath(self) -> str:
        """Local name of the DNA FASTA file."""
        return self._local_fasta_filepath(self.remote_dna_fasta_filename)

    @property
    def local_ncrna_fasta_filepath(self) -> str:
        """Local name of the ncRNA FASTA file."""
        return self._local_fasta_filepath(self.remote_ncrna_fasta_filename)

    @property
    def local_pep_fasta_filepath(self) -> str:
        """Local name of the peptide FASTA file."""
        return self._local_fasta_filepath(self.remote_pep_fasta_filename)

    def _local_fasta_filepath(self, remote_fasta_filename: str) -> str:
        return os.path.join(self.release_cache_dir, remote_fasta_filename)

    def download_cdna_fasta(self):
        """Download the cDNA FASTA file from the FTP server."""
        self._ftp_download(
            ENSEMBL_FTP_SERVER,
            self.remote_cdna_fasta_subdir,
            self.remote_cdna_fasta_filename,
            self.local_cdna_fasta_filepath,
        )

    def download_dna_fasta(self):
        """Download the DNA FASTA file from the FTP server."""
        self._ftp_download(
            ENSEMBL_FTP_SERVER,
            self.remote_dna_fasta_subdir,
            self.remote_dna_fasta_filename,
            self.local_dna_fasta_filepath,
        )

    def download_ncrna_fasta(self):
        """Download the ncRNA FASTA file from the FTP server."""
        self._ftp_download(
            ENSEMBL_FTP_SERVER,
            self.remote_ncrna_fasta_subdir,
            self.remote_ncrna_fasta_filename,
            self.local_ncrna_fasta_filepath,
        )

    def download_pep_fasta(self):
        """Download the peptide FASTA file from the FTP server."""
        self._ftp_download(
            ENSEMBL_FTP_SERVER,
            self.remote_pep_fasta_subdir,
            self.remote_pep_fasta_filename,
            self.local_pep_fasta_filepath,
        )

    @property
    def local_cdna_index_filepath(self) -> str:
        """Local name of the cDNA FASTA file."""
        return self._fasta_index_path(self.local_cdna_fasta_filepath)

    @property
    def local_dna_index_filepath(self) -> str:
        """Local name of the DNA FASTA file."""
        return self._fasta_index_path(self.local_dna_fasta_filepath)

    @property
    def local_ncrna_index_filepath(self) -> str:
        """Local name of the ncRNA FASTA file."""
        return self._fasta_index_path(self.local_ncrna_fasta_filepath)

    @property
    def local_pep_index_filepath(self) -> str:
        """Local name of the peptide FASTA file."""
        return self._fasta_index_path(self.local_pep_fasta_filepath)

    def _fasta_index_path(self, local_fasta_filename: str) -> str:
        """Return the path to the FASTA index file."""
        return local_fasta_filename + ".fai"

    def index_cdna_fasta(self):
        """(Re)build the index file for the cDNA Fasta."""
        return self._index_fasta(self.local_cdna_fasta_filepath)

    def index_dna_fasta(self):
        """(Re)build the index file for the DNA Fasta."""
        return self._index_fasta(self.local_dna_fasta_filepath)

    def index_ncrna_fasta(self):
        """(Re)build the index file for the cDNA Fasta."""
        return self._index_fasta(self.local_ncrna_fasta_filepath)

    def index_pep_fasta(self):
        """(Re)build the index file for the peptide Fasta."""
        return self._index_fasta(self.local_pep_fasta_filepath)

    def _index_fasta(self, local_fasta_filename: str):
        """(Re)build the index file for the FASTA file."""
        logger.info(f"Indexing {local_fasta_filename} (this may take some time)...")
        _ = Fasta(
            local_fasta_filename,
            key_function=strip_version,
            as_raw=True,
            sequence_always_upper=True,
            build_index=True,
            rebuild=True,
        )

    def load_cdna_fasta(self):
        """Load and return the cDNA Fasta."""
        return self._load_fasta(self.local_cdna_fasta_filepath)

    def load_dna_fasta(self):
        """Load and return the DNA Fasta."""
        return self._load_fasta(self.local_dna_fasta_filepath)

    def load_ncrna_fasta(self):
        """Load and return the cDNA Fasta."""
        return self._load_fasta(self.local_ncrna_fasta_filepath)

    def load_pep_fasta(self):
        """Load and return the peptide Fasta."""
        return self._load_fasta(self.local_pep_fasta_filepath)

    def _load_fasta(self, local_fasta_filename: str) -> Fasta:
        """Return a `pyfaidx.Fasta` object for the Ensembl genome."""
        return Fasta(
            local_fasta_filename,
            key_function=strip_version,
            as_raw=True,
            sequence_always_upper=True,
            build_index=False,
            rebuild=False,
        )

    # ---------------------------------------------------------------------------------------------
    # GTF file
    # ---------------------------------------------------------------------------------------------
    @property
    def remote_gtf_subdir(self) -> str:
        """Path to the directory on the FTP server that contains the GTF file."""
        return GTF_SUBDIR_TEMPLATE.format(release=self.release, species=self.species)

    @property
    def remote_gtf_filename(self) -> str:
        """Name of the GTF file on the FTP server."""
        return GTF_FILENAME_TEMPLATE.format(
            species=self.species.capitalize(), reference=self.reference, release=self.release
        )

    @property
    def local_gtf_filepath(self) -> str:
        """Local name of the GTF file."""
        return os.path.join(self.release_cache_dir, self.remote_gtf_filename)

    @property
    def local_gtf_cache_filepath(self) -> str:
        """Local name of the cache file."""
        return os.path.join(
            self.release_cache_dir, self.remote_gtf_filename.replace(".gtf.gz", ".pickle")
        )

    def download_gtf(self):
        """Download the GTF file from the FTP server."""
        self._ftp_download(
            ENSEMBL_FTP_SERVER,
            self.remote_gtf_subdir,
            self.remote_gtf_filename,
            self.local_gtf_filepath,
        )

    def normalize_gtf(self) -> pandas.DataFrame:
        """Normalize the data in the GTF to what's required by this package."""
        df = read_gtf(self.local_gtf_filepath)

        # # DEBUG
        # df = df[
        #     df["gene_name"].isin(
        #         [
        #             "AGBL1",
        #             "ALDOA",
        #             "BRCA2",
        #             "CHEK2",
        #             "KIT",
        #             "KMT2A",
        #             "MEN1",
        #             "MLL",
        #             "POU5F1",
        #             "PTEN",
        #             "SMARCA4",
        #             "SSX4",
        #             "TERT",
        #             "TCF19",
        #             "WNK1",
        #             "VPS53",
        #         ]
        #     )
        # ]
        # # /DEBUG

        logger.debug("_normalize_gtf_columns")
        df = self._normalize_gtf_columns(df)
        logger.debug("_add_missing_transcripts")
        df = self._add_missing_transcripts(df)
        logger.debug("_calculate_exon_offset_from_transcript")
        df = self._calculate_exon_offset_from_transcript(df)
        logger.debug("_calculate_cds_offset_from_transcript")
        df = self._calculate_cds_offset_from_transcript(df)
        logger.debug("_calculate_cds_offset_from_cdna")
        df = self._calculate_cds_offset_from_cdna(df)
        logger.debug("_assign_protein_id_to_stop_codon")
        df = self._assign_protein_id_to_stop_codon(df)

        # df.to_csv(self.remote_gtf_filename, sep="\t", index=False)  # DEBUG

        return df

    def _normalize_gtf_columns(self, df: pandas.DataFrame) -> pandas.DataFrame:
        """Rename columns, drop unused data, and set data types for each column."""
        # Rename column(s)
        df = df.rename(columns=GTF_COLUMN_RENAME)
        # Drop unused columns
        df = df.drop(df.columns.difference(GTF_KEEP_COLUMNS), axis=1)
        # Drop unused feature types
        df = df[df.feature.isin(GTF_KEEP_FEATURES)]
        # Coerce the non-null values in the 'exon_number' to integers
        df["exon_number"] = df["exon_number"].astype(int, errors="ignore")

        return df

    def _add_missing_transcripts(self, df: pandas.DataFrame) -> pandas.DataFrame:
        """Infer transcripts position(s) from other features, if not already defined."""
        transcript_rows = []

        def missing_transcript(group: pandas.DataFrame) -> pandas.Series:
            return group[group.feature == "transcript"].empty

        # add in missing transcripts
        for transcript_id, group in df.groupby("transcript_id"):
            if not transcript_id:
                continue

            if missing_transcript(group):
                features = group.sort_values(["start", "end"])
                first = features.iloc[0]
                last = features.iloc[-1]
                logger.debug(f"Adding transcript for '{transcript_id}'")
                transcript_rows.append(
                    {
                        "contig_id": first.contig_id,
                        "source": "",
                        "feature": "transcript",
                        "start": first.start,
                        "end": last.end,
                        "score": ".",
                        "strand": first.strand,
                        "frame": ".",
                        "gene_id": first.gene_id,
                        "gene_name": first.gene_name,
                        "transcript_id": transcript_id,
                        "transcript_name": first.transcript_name,
                    }
                )

        new_transcript_df = pandas.DataFrame(transcript_rows)
        df = pandas.concat([df, new_transcript_df], ignore_index=True)

        return df

    def _calculate_exon_offset_from_transcript(self, df: pandas.DataFrame) -> pandas.DataFrame:
        """Calculate the position of each exon, relative to the start of the transcript."""

        def get_transcript(group: pandas.DataFrame) -> pandas.Series:
            subdf = group[group.feature == "transcript"]
            assert len(subdf) == 1, transcript_id
            return subdf.iloc[0]

        result = {}

        for transcript_id, group in df.groupby("transcript_id"):
            if not transcript_id:
                continue

            offset = 0
            transcript = get_transcript(group)
            ascending = transcript.strand == "+"
            exon_df = group[group.feature == "exon"]
            for _, exon in exon_df.sort_values(["start", "end"], ascending=ascending).iterrows():
                # if this is the first exon/CDS/etc start the offset relative to the RNA
                if offset == 0:
                    if ascending:
                        offset = exon.start - transcript.start
                    else:
                        offset = transcript.end - exon.end
                key = (transcript_id, exon.feature, exon.start, exon.end, exon.strand)
                length = exon.end - exon.start + 1
                transcript_start = offset + 1
                transcript_end = length + offset
                result[key] = (transcript_start, transcript_end)
                offset += length

        # there's probably a better way to do this
        for index, feature in df.iterrows():
            key = (
                feature.transcript_id,
                feature.feature,
                feature.start,
                feature.end,
                feature.strand,
            )
            if (value := result.get(key)) is not None:  # type: ignore
                transcript_start, transcript_end = value
                df.loc[index, "transcript_start"] = transcript_start
                df.loc[index, "transcript_end"] = transcript_end

        # convert the offset values to integers
        df["transcript_start"] = df["transcript_start"].astype("Int64")
        df["transcript_end"] = df["transcript_end"].astype("Int64")

        return df

    def _calculate_cds_offset_from_transcript(self, df: pandas.DataFrame) -> pandas.DataFrame:
        """Calculate the position of each CDS, relative to the start of the transcript."""

        def get_exon(
            group: pandas.DataFrame, transcript_id: str, start: int, end: int
        ) -> pandas.Series:
            mask = (group.start <= start) & (group.end >= end) & (group.feature == "exon")
            subdf = group[mask]
            assert len(subdf) == 1, (transcript_id, start, end)
            return subdf.iloc[0]

        result = {}

        for transcript_id, group in df.groupby("transcript_id"):
            cds_df = group[group.feature.isin(["CDS", "stop_codon"])]
            if cds_df.empty:
                continue

            ascending = cds_df.iloc[0].strand == "+"
            for _, cds in cds_df.sort_values(["start", "end"], ascending=ascending).iterrows():
                exon = get_exon(group, transcript_id, cds.start, cds.end)
                if ascending:
                    offset = cds.start - exon.start + exon.transcript_start
                else:
                    offset = exon.end - cds.end + exon.transcript_start
                key = (transcript_id, cds.feature, cds.start, cds.end, cds.strand)
                length = cds.end - cds.start
                transcript_start = offset
                transcript_end = length + offset
                result[key] = (transcript_start, transcript_end)
                offset += length

        # there's probably a better way to do this
        for index, feature in df.iterrows():
            key = (
                feature.transcript_id,
                feature.feature,
                feature.start,
                feature.end,
                feature.strand,
            )
            if (value := result.get(key)) is not None:  # type: ignore
                transcript_start, transcript_end = value
                df.loc[index, "transcript_start"] = transcript_start
                df.loc[index, "transcript_end"] = transcript_end

        # convert the offset values to integers
        df["transcript_start"] = df["transcript_start"].astype("Int64")
        df["transcript_end"] = df["transcript_end"].astype("Int64")

        return df

    def _calculate_cds_offset_from_cdna(self, df: pandas.DataFrame) -> pandas.DataFrame:
        """Calculate the position of each CDS, relative to the start of the cDNA."""

        def infer_cdna(subdf: pandas.DataFrame):
            if len(subdf) > 0:
                subdf_sorted = subdf.sort_values(["start", "end"], ascending=True)
                ascending = subdf.iloc[0].strand == "+"
                cdna_start = subdf_sorted.iloc[0].start
                cdna_end = subdf_sorted.iloc[-1].end
                return cdna_start, cdna_end, ascending
            else:
                return None

        result = {}

        for transcript_id, group in df.groupby("transcript_id"):
            if not transcript_id:
                continue

            cds_df = group[group.feature.isin(["CDS", "stop_codon"])]
            if cdna := infer_cdna(cds_df):
                cdna_start_pos, cdna_end_pos, ascending = cdna
            else:
                continue

            offset = 0
            for _, cds in cds_df.sort_values(["start", "end"], ascending=ascending).iterrows():
                # if this is the first exon/CDS/etc start the offset relative to the RNA
                if offset == 0:
                    if ascending:
                        offset = cds.start - cdna_start_pos
                    else:
                        offset = cdna_end_pos - cds.end
                key = (transcript_id, cds.feature, cds.start, cds.end, cds.strand)
                length = cds.end - cds.start + 1
                cdna_start = offset + 1
                cdna_end = length + offset
                result[key] = (cdna_start, cdna_end)
                offset += length

        # there's probably a better way to do this
        for index, feature in df.iterrows():
            key = (
                feature.transcript_id,
                feature.feature,
                feature.start,
                feature.end,
                feature.strand,
            )
            if (value := result.get(key)) is not None:  # type: ignore
                cdna_start, cdna_end = value
                df.loc[index, "cdna_start"] = cdna_start
                df.loc[index, "cdna_end"] = cdna_end

        # convert the offset values to integers
        df["cdna_start"] = df["cdna_start"].astype("Int64")
        df["cdna_end"] = df["cdna_end"].astype("Int64")

        return df

    def _assign_protein_id_to_stop_codon(self, df: pandas.DataFrame) -> pandas.DataFrame:
        """Add the protein ID as info for each stop codon, if not already defined."""

        def select(group: pandas.DataFrame) -> pandas.Series:
            subdf = group[group.feature == "CDS"]
            return subdf.iloc[0] if len(subdf) > 0 else None

        for _, group in df.groupby("transcript_id"):
            if (cds := select(group)) is None:
                continue

            stop_codon_df = group[group.feature == "stop_codon"]
            for index, _ in stop_codon_df.iterrows():
                if not df.loc[index, "protein_id"]:
                    df.loc[index, "protein_id"] = cds.protein_id

        return df

    def cache_gtf(self, df: pandas.DataFrame) -> pandas.DataFrame:
        """Convert the Ensembl GTF to a pandas DataFrame and cache it."""
        logger.debug(f"Converting {self.local_gtf_filepath} to {self.local_gtf_cache_filepath}")
        df.to_pickle(self.local_gtf_cache_filepath)
        # logger.debug(f"Removing {self.local_gtf_filepath}")
        # os.remove(self.local_gtf_filepath)  # DEBUG

    def load_gtf(self) -> pandas.DataFrame:
        """Return the cached Ensembl GTF data."""
        return pandas.read_pickle(self.local_gtf_cache_filepath)

    # ---------------------------------------------------------------------------------------------
    # generic functions
    # ---------------------------------------------------------------------------------------------
    def _ftp_download(self, server: str, subdir: str, remote_file: str, local_file: str):
        """Download a file from an FTP server."""
        try:
            ftp = FTP(server)
            logger.debug(f"Connecting to {server}")
            ftp.login()
            url = "ftp://" + server + "/" + subdir + "/" + remote_file
            logger.debug(f"Downloading {url} to {local_file}")
            ftp.cwd(subdir)
            self.make_release_cache_dir()
            with open(local_file, "wb") as fh:
                ftp.retrbinary(f"RETR {remote_file}", fh.write)
            logger.debug("Download successful")
            ftp.quit()
        except Exception as exc:
            logger.exception(f"Download failed: {exc}")
