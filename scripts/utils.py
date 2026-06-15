"""
Utility functions for the GROMACS simulation workflow.
"""

import re
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
        if re.match(r"#\s+gmx\s+distance\s+-f", line):
            try:
                select = line.split("-select '")[1].split("' -")[0]
            except IndexError as e:
                raise IndexError("Failed to parse select string from line: " + line) from e

            columns = select.split("; ")
        if re.match(r"@\s+xaxis\s+label", line):
            x_unit = line.split("label \"")[1].split("\"")[0]
        if re.match(r"@\s+yaxis\s+label", line):
            y_unit = line.split("label \"")[1].split("\"")[0]

    if x_unit is None:
        raise ValueError("Could not find x-axis unit in XVG file")
    if y_unit is None:
        raise ValueError("Could not find y-axis unit in XVG file")
    if columns is None:
        raise ValueError("Could not find columns in XVG file")

    columns = [x_unit] + [f"{column} {y_unit}" for column in columns]

    data = "\n".join(
        line for line in data.splitlines()
        if line.strip() and not line.startswith(("#", "@", "&"))
    )

    print(columns)
    print(data[:200])
    df = pd.read_csv(StringIO(data), sep=r"\s+", names=columns)

    return df