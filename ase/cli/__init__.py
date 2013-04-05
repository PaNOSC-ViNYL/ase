import sys
import traceback

from ase.cli.cli import run


def main(default_calculator={'name': 'emt'}):
    try:
        run(default_calculator=default_calculator)
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
