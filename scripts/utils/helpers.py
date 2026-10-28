import logging
import os
import re
import sys
from pathlib import Path

def gen_flags(args):
    """ Generate command-line flags from argument dictionary. """
    flags = []

    for name, value in args.items():
        if value is not None:
            if isinstance(value, bool) and not value:
                continue
            flag = "--" + name
            flags.append(flag)
            if type(value) != bool:
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



# Environment variable that makes the console show DEBUG messages too (set by main.py --verbose).
VERBOSE_ENV = "MDPIPE_VERBOSE"


def console_verbose() -> bool:
    """True if DEBUG messages should also be printed on the console."""
    return os.environ.get(VERBOSE_ENV) == "1"


def setup_logger(name: str, logfile: Path, level=logging.DEBUG):
    """Logger writing all messages to ``logfile`` and INFO and above to stdout.

    With the environment variable MDPIPE_VERBOSE=1 the console also shows DEBUG messages.
    The first call for a logger name sets it up; later calls return it unchanged.
    """
    logger = logging.getLogger(name)

    if logger.handlers:
        return logger

    logger.setLevel(level)

    # The log file receives every message; the console only INFO and above unless verbose.
    fh = logging.FileHandler(str(logfile), mode="a")
    fh.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG if console_verbose() else logging.INFO)

    formatter = logging.Formatter(
        "%(levelname)s - %(message)s"
    )

    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger
