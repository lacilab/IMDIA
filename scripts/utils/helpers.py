import logging
import os
import re


# Generate command-line flags from argument dictionary
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

def setup_logger(name, logfile, level=logging.DEBUG):

    logger = logging.getLogger(name)

    if logger.handlers:
        return logger

    os.makedirs(os.path.dirname(logfile), exist_ok=True)

    logger.setLevel(level)

    fh = logging.FileHandler(logfile, mode="a")
    ch = logging.StreamHandler()

    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s"
    )

    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger
