import sys

import numpy as np
import pandas as pd

from .constants import CONTIG_ID

GTF_COLUMN_RENAME = {"seqname": CONTIG_ID}
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
    df = df.sort_values(["start", "end"])
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
    assert sorted(df.columns) == sorted(GTF_KEEP_COLUMNS)
    # Coerce the non-null values in the 'exon_number' to float
    df["exon_number"] = df["exon_number"].astype(float, errors="raise")

    return df


def df_infer_transcripts(df: pd.DataFrame) -> pd.DataFrame:
    """Infer transcripts position(s) from other features, if not already defined."""
    transcript_rows = []

    for transcript_id, group in df.groupby("transcript_id"):
        if group[group.feature == "transcript"].empty:
            first = group.iloc[0]
            last = group.iloc[-1]
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
    df = df.sort_values(["start", "end"])

    return df


def df_infer_cdna(df: pd.DataFrame) -> pd.DataFrame:
    """Infer cDNA position(s) from other features, if not already defined."""
    cdna_rows = []

    for transcript_id, group in df.groupby("transcript_id"):
        if group[group.feature == "cdna"].empty:
            cds_df = group[group.feature == "CDS"]
            if not cds_df.empty:
                first = cds_df.iloc[0]
                last = cds_df.iloc[-1]
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
    df = df.sort_values(["start", "end"])

    return df


def df_infer_missing_genes(df: pd.DataFrame) -> pd.DataFrame:
    """Infer gene position(s) from other features, if not already defined."""
    gene_rows = []

    for gene_id, group in df.groupby("gene_id"):
        if group[group.feature == "gene"].empty:
            first = group.iloc[0]
            last = group.iloc[-1]
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
    df = df.sort_values(["start", "end"])

    return df


def df_exon_offset_transcript(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate the position of each exon, relative to the start of the transcript."""
    for transcript_id, group in df.groupby("transcript_id"):
        transcript_df = group[group.feature == "transcript"]
        if transcript_df.empty:
            continue

        offset = 0
        transcript_index = transcript_df.index[0]
        transcript = transcript_df.iloc[0]
        ascending = transcript.strand == "+"
        exon_df = group[group.feature == "exon"]
        if exon_df.empty:
            continue

        for exon in exon_df.sort_values(["start", "end"], ascending=ascending).itertuples():
            # if this is the first exon/CDS/etc start the offset relative to the RNA
            if offset == 0:
                if ascending:
                    offset = exon.start - transcript.start
                else:
                    offset = transcript.end - exon.end

            length = exon.end - exon.start + 1
            transcript_start = offset + 1
            transcript_end = length + offset
            offset += length

            # add the new offsets to the exon
            df.at[exon.Index, "transcript_start"] = transcript_start
            df.at[exon.Index, "transcript_end"] = transcript_end

        # add the new offsets to the transcript
        df.at[transcript_index, "transcript_start"] = 1
        df.at[transcript_index, "transcript_end"] = offset

    # convert the offset values to integers
    df["transcript_start"] = df["transcript_start"].astype("Int64")
    df["transcript_end"] = df["transcript_end"].astype("Int64")

    return df


def df_cds_offset_transcript(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate the position of each CDS, relative to the start of the transcript."""
    for transcript_id, group in df.groupby("transcript_id"):
        cds_df = group[group.feature.isin(["CDS", "stop_codon"])]
        if cds_df.empty:
            continue

        ascending = cds_df.iloc[0].strand == "+"
        for cds in cds_df.sort_values(["start", "end"], ascending=ascending).itertuples():
            exon_mask = (
                (group.start <= cds.start) & (group.end >= cds.end) & (group.feature == "exon")
            )
            exon_df = group[exon_mask]
            if exon_df.empty:
                continue

            exon = exon_df.iloc[0]
            if ascending:
                offset = cds.start - exon.start + exon.transcript_start
            else:
                offset = exon.end - cds.end + exon.transcript_start

            length = cds.end - cds.start
            transcript_start = offset
            transcript_end = length + offset
            offset += length

            # add the new offsets and exon ID to the CDS
            df.at[cds.Index, "exon_id"] = exon.exon_id
            df.at[cds.Index, "transcript_start"] = transcript_start
            df.at[cds.Index, "transcript_end"] = transcript_end

    # convert the offset values to integers
    df["transcript_start"] = df["transcript_start"].astype("Int64")
    df["transcript_end"] = df["transcript_end"].astype("Int64")

    return df


def df_cds_offset_cdna(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate the position of each CDS, relative to the start of the cDNA."""
    for transcript_id, group in df.groupby("transcript_id"):
        cdna_df = group[group.feature == "cdna"]
        if cdna_df.empty:
            continue

        cdna_index = cdna_df.index[0]
        cdna = cdna_df.iloc[0]
        ascending = cdna.strand == "+"

        cds_df = group[group.feature.isin(["CDS", "stop_codon"])]
        if cds_df.empty:
            continue

        offset = 0
        cds_df = cds_df.sort_values(["start", "end"], ascending=ascending)
        for cds in cds_df.itertuples():
            # if this is the first exon/CDS/etc start the offset relative to the RNA
            if offset == 0:
                if ascending:
                    offset = cds.start - cdna.start
                else:
                    offset = cdna.end - cds.end

            length = cds.end - cds.start + 1
            start_offset = offset + 1
            end_offset = length + offset
            offset += length

            # add the new offsets to the CDS
            df.at[cds.Index, "cdna_start"] = start_offset
            df.at[cds.Index, "cdna_end"] = end_offset

        # add the new offsets to the cDNA
        df.at[cdna_index, "cdna_start"] = 1
        df.at[cdna_index, "cdna_end"] = offset

    # convert the offset values to integers
    df["cdna_start"] = df["cdna_start"].astype("Int64")
    df["cdna_end"] = df["cdna_end"].astype("Int64")

    return df


def df_set_protein_id(df: pd.DataFrame) -> pd.DataFrame:
    """Add the protein ID as info for each cDNA and stop codon, if not already defined."""
    for _, group in df.groupby("transcript_id"):
        # Get the first non-NA protein ID matching the transcript ID
        # A protein coding transcript should only have 1 matching protein ID
        if pidx := group["protein_id"].first_valid_index():
            protein_id = group["protein_id"].loc[pidx]
        else:
            continue

        subdf = group[group.feature.isin(["cdna", "stop_codon"])]
        for index, _ in subdf.iterrows():
            if pd.isna(df.loc[index, "protein_id"]):
                print(f"Adding protein ID '{protein_id}' to row {index}", file=sys.stderr)
                df.loc[index, "protein_id"] = protein_id

    return df
