# libcluster.py
import socket
import re

def detect_cluster():
    hostname = socket.gethostname().lower()

    if re.match(r"login0[1-4](?:\.|$)", hostname):
        return "kisti"
    elif hostname.startswith("tgm-master"):
        return "pt"
    else:
        return "unknown"
