"""
Utility functions for the GROMACS simulation workflow.
"""

import pandas as pd
from io import StringIO
from pathlib import Path

def read_xvg(file_path):
    """
    Read an XVG file created with a -select option and return a pandas DataFrame.
    """
    x_unit = None
    y_unit = None
    columns = None

    path = Path(file_path)

    data = path.read_text(encoding="utf-8")

    for line in data.splitlines():
        if line.startswith("# gmx distance"):
            select = line.split("-select '")[1].split("' -")[0]
            columns = select.split("; ")
        if line.startswith("@    xaxis  label"):
            x_unit = line.split("label \"")[1].split("\"")[0]
        if line.startswith("@    yaxis  label"):
            y_unit = line.split("label \"")[1].split("\"")[0]

    assert (x_unit is not None and
            y_unit is not None and
            columns is not None), "Failed to parse XVG file"

    columns = [x_unit] + [f"{column} {y_unit}" for column in columns]

    data = "\n".join(
        line for line in data.splitlines()
        if line.strip() and not line.startswith(("#", "@", "&"))
    )

    df = pd.read_csv(StringIO(data), sep=r"\s+", names=columns)

    return df