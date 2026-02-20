import socket

def detect_cluster():
    hostname = socket.gethostname()

    if "pt" in hostname:
        return "pt"
    if "kisti" in hostname:
        return "kisti"

    return "local"
