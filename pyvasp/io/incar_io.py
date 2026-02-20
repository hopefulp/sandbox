def modify_incar(path, changes, remove=None):

    with open(path, "r") as f:
        lines = f.readlines()

    new_lines = []

    for line in lines:
        key = line.split("=")[0].strip()

        if remove and key in remove:
            continue

        if key in changes:
            new_lines.append(f"{key} = {changes[key]}\n")
            del changes[key]
        else:
            new_lines.append(line)

    for k, v in changes.items():
        new_lines.append(f"{k} = {v}\n")

    with open("INCAR.new", "w") as f:
        f.writelines(new_lines)
