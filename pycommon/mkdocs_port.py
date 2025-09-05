#!/home/joonho/anaconda3/bin/python
import socket
import subprocess
import webbrowser

# Find an open port
sock = socket.socket()
sock.bind(('', 0))  # bind to a random available port
port = sock.getsockname()[1]
webbrowser.open(f"http://127.0.0.1:{port}")
sock.close()

# Run mkdocs serve with that port
subprocess.run(["mkdocs", "serve", "--dev-addr", f"127.0.0.1:{port}"])

