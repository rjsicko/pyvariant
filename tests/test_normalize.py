import numpy as np
import pandas as pd

from ensembl_map.normalize import df_set_protein_id


def test_df_set_protein_id():
    columns = ["feature", "transcript_id", "protein_id"]
    data = [
        ["CDS", "ENST00000123456", "ENSP00000654321"],
        ["CDS", "ENST00000412167", "ENSP00000390987"],
        ["exon", "ENST00000999999", np.nan],
        ["stop_codon", "ENST00000412167", np.nan],
        ["cdna", "ENST00000412167", np.nan],
    ]
    in_df = pd.DataFrame(data, columns=columns)
    out_df = df_set_protein_id(in_df)
    assert out_df["protein_id"].to_list() == [
        "ENSP00000654321",
        "ENSP00000390987",
        np.nan,
        "ENSP00000390987",
        "ENSP00000390987",
    ]
