import sys
import traceback

from ase.cli.cli import run


def main(hook=None):
    try:
        run(hook=hook)
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception:
        traceback.print_exc()
        sys.stderr.write("""
An exception occurred!  Please report the issue to
ase-developer@listserv.fysik.dtu.dk - thanks!  Please also report this
if it was a user error, so that a better error message can be provided
next time.""")
        raise
