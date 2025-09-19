'''
data-file related functions
'''
import pandas as pd
from libstr import isnumber, find_sep

def load_table(fin, icx, jcy, sep):
    """
    data file has table
    Extract tabular data from file with x in i-th ys in js using pandas.
    
    Parameters:
        filepath (str): Path to the file.
        icx (int): Index of column for x-axis with index starting from 0
        jcy (int): Index of column for y-axis with index starting from 0
        sep         auto-detect
                (old) sep (str, optional): Delimiter, e.g., ',' for CSV, '\t' for TSV.
                             If None, pandas will try to guess.
    """

    # First read one row to check file format and header 
    with open(fin, "r") as f:
        first_line = f.readline().strip()
        #second_line = f.readline().strip()
        sep = find_sep(first_line)
        ### default is delim_whitespace
    first_list = first_line.split(sep)
    #print(f"1st line {first_line} and list {first_list}")

    # Decide if header exists (if the j-th element is not numeric)
    
    #print(f"index jcy[0] {jcy[0]} first_list {first_list}")
    has_header = not isnumber(first_list[jcy[0]])

    # Read file accordingly
    if sep is None:
        df = pd.read_csv(fin, delim_whitespace=True, header=0 if has_header else None)
    else:
        df = pd.read_csv(fin, sep=sep, header=0 if has_header else None)

    x = df.iloc[:, icx]
    xlabel = df.columns[icx] if has_header else f"Column {icx}"

    ys, yls = [], []

    for j in jcy:
        ys.append(df.iloc[:, j])
        yls.append(df.columns[j] if has_header else f"Column {j}")

    return x, xlabel, ys, yls