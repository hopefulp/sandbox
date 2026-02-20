# libcluster.py
import socket
import re

def detect_cluster():
    hostname = socket.gethostname()

    if hostname.startswith("login0"):
        return "kisti"
    elif hostname.startswith("tgm-master"):
        return "pt"
    else:
        return "unknown"
