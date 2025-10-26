import os

def fname_alpha_contract(absdir):
    """
    Convert path like H2p_00.70/3.alpha/rc4.80
    → jobname like d0070rc480
    """
    base = os.path.basename(absdir)          # e.g. "rc4.80"
    parent = os.path.basename(os.path.dirname(absdir))  # e.g. "3.alpha"
    grandparent = os.path.basename(os.path.dirname(os.path.dirname(absdir)))  # e.g. "H2p_00.70"

    # from "H2p_00.70" → "0070"
    dist = grandparent.split("_")[-1]  # "00.70"
    dist = dist.replace(".", "").zfill(4)  # "0070"

    # from "rc4.80" → "rc480"
    rc = base.replace(".", "")  # "rc480"

    return f"d{dist}{rc}"
