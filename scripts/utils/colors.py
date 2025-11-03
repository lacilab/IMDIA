import sys

# Color definitions
RED = "\033[0;31m"
GREEN = "\033[0;32m"
YELLOW = "\033[0;33m"
NC = "\033[0m"  # No color / reset

# Helper print functions
def ok(msg):
    print(f"{GREEN}[OK]{NC} {msg}")

def warn(msg):
    print(f"{YELLOW}[WARNING]{NC} {msg}")

def error(msg):
    print(f"{RED}[ERROR]{NC} {msg}", file=sys.stderr)
