# Generate command-line flags from argument dictionary
import re


def gen_flags(args):
    flags = []

    for name, value in args.items():
        if value != False and value is not None:
            flag = "--" + name
            flags.append(flag)
            if value != True:
                flags.append(str(value))

    return flags

def count_lines(file, *names):
    counts = {name: 0 for name in names}

    patterns = {
        name: re.compile(rf'(?<!\S){re.escape(name)}(?!\S)')
        for name in names
    }

    with open(file) as f:
        for line in f:
            for name, pattern in patterns.items():
                counts[name] += len(pattern.findall(line))

    return [counts[name] for name in names]