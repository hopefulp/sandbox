# ============================================================
# Generic utilities
# ============================================================

import os

def walk_dirs( root, dir_work, *, follow_symlinks=False, exclude_dirs=None, **work_kwargs ):
    """
    Generic recursive directory walker.
    Calls dir_work(path, **work_kwargs) in every directory.
        dir_clean
        dir_touch
        etc
    """
    if not root:
        root = os.getcwd()

    if exclude_dirs is None:
        exclude_dirs = set()

    root = os.path.abspath(root)

    print(f"#### enter {root}")

    # DO WORK IN THIS DIRECTORY
    dir_work(root, **work_kwargs)

    # RECURSE
    for name in os.listdir(root):
        path = os.path.join(root, name)

        if not os.path.isdir(path):
            continue
        if not follow_symlinks and os.path.islink(path):
            continue
        if name in exclude_dirs:
            continue

        walk_dirs(
            path,
            dir_work,
            follow_symlinks=follow_symlinks,
            exclude_dirs=exclude_dirs,
            **work_kwargs,
        )

    print(f"#### exit  {root}")
