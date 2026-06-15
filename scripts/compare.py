from functools import reduce
from glob import glob
from pathlib import Path

import pandas as pd

from scripts import utils

def combine_xvg(file_glob : str, x_column : str, y_columns : list[str] = None) -> pd.DataFrame:
    """
    Combine a list of XVG files into a single DataFrame. Taking specified columns from each file.
    """
    dataframes = []
    for file in glob(file_glob):
        df = utils.read_xvg(file)

        if y_columns is not None:
            df = df[[x_column] + y_columns]

        # Suffix every column except the shared x column with the file name
        # so columns coming from different files stay distinct after merging.
        suffix = Path(file).stem
        df = df.rename(columns={
            col: f"{col}_{suffix}"
            for col in df.columns
            if col != x_column
        })

        dataframes.append(df)

    return reduce(
        lambda left, right: pd.merge(left, right, on=x_column, how="outer"),
        dataframes,
    )

def filter_columns(df : pd.DataFrame, contains : list[str]) -> pd.DataFrame:
    """
    Filter columns of a DataFrame by a list of strings.
    """
    return df.filter(regex=f"({'|'.join(contains)})")