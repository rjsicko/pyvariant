import sys

import numpy as np
import pandas as pd

GTF_COLUMN_RENAME = {"seqname": "contig_id"}
GTF_KEEP_FEATURES = ["CDS", "exon", "gene", "stop_codon", "transcript"]
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


def normalize_df(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize a pandas DataFrame and infer missing data."""
    df = df.replace("", np.nan)
    df = df_normalize_cols(df)
    df = df_infer_transcripts(df)
    df = df_infer_cdna(df)
    df = df_infer_missing_genes(df)
    df = df_exon_offset_transcript(df)
    df = df_cds_offset_transcript(df)
    df = df_cds_offset_cdna(df)
    df = df_set_protein_id(df)
    df = df.replace("", np.nan)

    return df


def df_normalize_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Rename columns, drop unused data, and set data types for each column."""
    # Rename column(s)
    df = df.rename(columns=GTF_COLUMN_RENAME)
    # Drop unused columns
    df = df.drop(df.columns.difference(GTF_KEEP_COLUMNS), axis=1)
    # Drop unused feature types
    df = df[df.feature.isin(GTF_KEEP_FEATURES)]
    # Assert that all the expected columns exist
    assert df.columns == GTF_KEEP_FEATURES
    # Coerce the non-null values in the 'exon_number' to integers
    df["exon_number"] = df["exon_number"].astype(int, errors="ignore")

    return df


def df_infer_transcripts(df: pd.DataFrame) -> pd.DataFrame:
    """Infer transcripts position(s) from other features, if not already defined."""
    transcript_rows = []

    def missing_transcript(group: pd.DataFrame) -> pd.Series:
        return group[group.feature == "transcript"].empty

    # add in missing transcripts
    for transcript_id, group in df.groupby("transcript_id"):
        if not transcript_id:
            continue

        if missing_transcript(group):
            features = group.sort_values(["start", "end"])
            first = features.iloc[0]
            last = features.iloc[-1]
            new_row = {
                "contig_id": first.contig_id,
                "feature": "transcript",
                "start": first.start,
                "end": last.end,
                "strand": first.strand,
                "gene_id": first.gene_id,
                "gene_name": first.gene_name,
                "transcript_id": transcript_id,
                "transcript_name": first.transcript_name,
            }
            transcript_rows.append(new_row)
            print(f"Inferred transcript {new_row}", file=sys.stderr)

    new_transcript_df = pd.DataFrame(transcript_rows)
    df = pd.concat([df, new_transcript_df], ignore_index=True)

    return df


def df_infer_cdna(df: pd.DataFrame) -> pd.DataFrame:
    """Infer cDNA position(s) from other features, if not already defined."""
    cdna_rows = []

    def missing_cdna(group: pd.DataFrame) -> pd.Series:
        return group[group.feature == "cdna"].empty

    # add in missing cDNA
    for transcript_id, group in df.groupby("transcript_id"):
        if not transcript_id:
            continue

        if missing_cdna(group):
            cds_df = group[group.feature == "CDS"]
            if not cds_df.empty:
                features = cds_df.sort_values(["start", "end"])
                first = features.iloc[0]
                last = features.iloc[-1]
                new_row = {
                    "contig_id": first.contig_id,
                    "feature": "cdna",
                    "start": first.start,
                    "end": last.end,
                    "strand": first.strand,
                    "gene_id": first.gene_id,
                    "gene_name": first.gene_name,
                    "transcript_id": transcript_id,
                    "transcript_name": first.transcript_name,
                }
                cdna_rows.append(new_row)
                print(f"Inferred cDNA {new_row}", file=sys.stderr)

    new_cdna_df = pd.DataFrame(cdna_rows)
    df = pd.concat([df, new_cdna_df], ignore_index=True)

    return df


def df_infer_missing_genes(df: pd.DataFrame) -> pd.DataFrame:
    """Infer gene position(s) from other features, if not already defined."""
    gene_rows = []

    def missing_gene(group: pd.DataFrame) -> pd.Series:
        return group[group.feature == "gene"].empty

    # add in missing transcripts
    for gene_id, group in df.groupby("gene_id"):
        if not gene_id:
            continue

        if missing_gene(group):
            features = group.sort_values(["start", "end"])
            first = features.iloc[0]
            last = features.iloc[-1]
            new_row = {
                "contig_id": first.contig_id,
                "feature": "gene",
                "start": first.start,
                "end": last.end,
                "strand": first.strand,
                "gene_id": first.gene_id,
                "gene_name": first.gene_name,
            }
            gene_rows.append(new_row)
            print(f"Inferred gene {new_row}", file=sys.stderr)

    new_gene_df = pd.DataFrame(gene_rows)
    df = pd.concat([df, new_gene_df], ignore_index=True)

    return df


def df_exon_offset_transcript(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate the position of each exon, relative to the start of the transcript."""

    def get_transcript(group: pd.DataFrame) -> pd.Series:
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

        # add in the transcript length to the transcript itself
        key = (
            transcript_id,
            transcript.feature,
            transcript.start,
            transcript.end,
            transcript.strand,
        )
        result[key] = (1, offset)

    # there's probably a better way to do this
    for index, feature in df.iterrows():
        key = (feature.transcript_id, feature.feature, feature.start, feature.end, feature.strand)
        if (value := result.get(key)) is not None:  # type: ignore
            transcript_start, transcript_end = value
            df.loc[index, "transcript_start"] = transcript_start
            df.loc[index, "transcript_end"] = transcript_end

    # convert the offset values to integers
    df["transcript_start"] = df["transcript_start"].astype("Int64")
    df["transcript_end"] = df["transcript_end"].astype("Int64")

    return df


def df_cds_offset_transcript(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate the position of each CDS, relative to the start of the transcript."""

    def get_exon(group: pd.DataFrame, transcript_id: str, start: int, end: int) -> pd.Series:
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
            result[key] = (transcript_start, transcript_end, exon.exon_id)
            offset += length

    # there's probably a better way to do this
    for index, feature in df.iterrows():
        key = (feature.transcript_id, feature.feature, feature.start, feature.end, feature.strand)
        if (value := result.get(key)) is not None:  # type: ignore
            transcript_start, transcript_end, exon_id = value
            df.loc[index, "exon_id"] = exon_id
            df.loc[index, "transcript_start"] = transcript_start
            df.loc[index, "transcript_end"] = transcript_end

    # convert the offset values to integers
    df["transcript_start"] = df["transcript_start"].astype("Int64")
    df["transcript_end"] = df["transcript_end"].astype("Int64")

    return df


def df_cds_offset_cdna(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate the position of each CDS, relative to the start of the cDNA."""
    result = {}

    for transcript_id, group in df.groupby("transcript_id"):
        if not transcript_id:
            continue

        cds_df = group[group.feature.isin(["CDS", "stop_codon"])]
        if cds_df.empty:
            continue

        cdna_df = group[group.feature == "cdna"]
        if cdna_df.empty:
            print(f"No cDNA found for transcript {transcript_id}", file=sys.stderr)
            continue
        else:
            if len(cdna_df) > 1:
                print(f"Multiple cDNA found for transcript {transcript_id}", file=sys.stderr)

            cdna = cdna_df.iloc[0]
            ascending = cdna.strand == "+"

        cds_df = group[group.feature.isin(["CDS", "stop_codon"])]
        if cds_df.empty:
            continue

        offset = 0
        cds_df = cds_df.sort_values(["start", "end"], ascending=ascending)
        for _, cds in cds_df.iterrows():
            # if this is the first exon/CDS/etc start the offset relative to the RNA
            if offset == 0:
                if ascending:
                    offset = cds.start - cdna.start
                else:
                    offset = cdna.end - cds.end
            key = (transcript_id, cds.feature, cds.start, cds.end, cds.strand)
            length = cds.end - cds.start + 1
            start_offset = offset + 1
            end_offset = length + offset
            result[key] = (start_offset, end_offset)
            offset += length

        # add in the cDNA length to the cDNA itself
        key = (transcript_id, cdna.feature, cdna.start, cdna.end, cdna.strand)
        result[key] = (1, offset)

    # there's probably a better way to do this
    for index, feature in df.iterrows():
        key = (feature.transcript_id, feature.feature, feature.start, feature.end, feature.strand)
        if (value := result.get(key)) is not None:  # type: ignore
            start_offset, end_offset = value
            df.loc[index, "cdna_start"] = start_offset
            df.loc[index, "cdna_end"] = end_offset

    # convert the offset values to integers
    df["cdna_start"] = df["cdna_start"].astype("Int64")
    df["cdna_end"] = df["cdna_end"].astype("Int64")

    return df


def df_set_protein_id(df: pd.DataFrame) -> pd.DataFrame:
    """Add the protein ID as info for each cDNA and stop codon, if not already defined."""

    def select(group: pd.DataFrame) -> pd.Series:
        subdf = group[group.feature == "CDS"]
        return subdf.iloc[0] if len(subdf) > 0 else None

    for _, group in df.groupby("transcript_id"):
        if (cds := select(group)) is None:
            continue

        subdf = group[group.feature.isin(["cdna", "stop_codon"])]
        # subdf.iloc[:, 'protein_id'].fillna(cds.protein_id, inplace=True)
        for index, _ in subdf.iterrows():
            if pd.isna(df.loc[index, "protein_id"]):
                print(f"Adding protein ID '{cds.protein_id}' to row {index}", file=sys.stderr)
                df.loc[index, "protein_id"] = cds.protein_id

    return df
