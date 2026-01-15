import os

def dir_touch(path):
    print(f"[TOUCH] {path}")

    for name in os.listdir(path):
        fpath = os.path.join(path, name)
        if os.path.isfile(fpath) and not os.path.islink(fpath):
            os.utime(fpath, None)
