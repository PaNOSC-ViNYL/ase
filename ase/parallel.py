import sys
import time
import atexit

def paropen(name, mode='r', buffering=0):
    """MPI-safe version of open function."""
    if rank > 0 and mode[0] != 'r':
        name = '/dev/null'
    return open(name, mode, buffering)


# Check for special MPI-enabled Python interpreters:
if '_gpaw' in sys.modules:
    # http://wiki.fysik.dtu.dk/gpaw
    from gpaw.mpi import world
    rank = world.rank
    size = world.size
    barrier = world.barrier
elif 'Scientific_mpi' in sys.modules:
    # 
    from Scientific.MPI import world
    rank = world.rank
    size = world.size
    barrier = world.barrier
else:
    # This is a standard Python interpreter:
    rank = 0
    size = 1
    def barrier():
        pass


def register_parallel_cleanup_function():
    """Call MPI_Abort if python crashes.

    This will terminate the processes on the other nodes."""
        
    if size == 1:
        return

    def cleanup(sys=sys, time=time, world=world):
        error = getattr(sys, 'last_type', None)
        if error:
            sys.stdout.flush()
            sys.stderr.write(('ASE CLEANUP (node %d): %s occurred.  ' +
                              'Calling MPI_Abort!\n') % (world.rank, error))
            sys.stderr.flush()
            # Give other nodes a moment to crash by themselves (perhaps
            # producing helpful error messages):
            time.sleep(3)
            world.abort(42)

    atexit.register(cleanup)
