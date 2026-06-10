import dask
# import dask.multiprocessing
from MDAnalysis import Universe
# from dask.distributed import Client
import numpy as np

def setup_dask(n_jobs: int):
    """
    Configure Dask scheduler for parallel execution.
    Uses multiprocessing to distribute tasks across n_jobs workers.

    Arguments:
        n_jobs : Number of parallel workers.
    """
    # client = Client(n_workers=n_jobs)
    dask.config.set(scheduler='processes')

@dask.delayed
def analyze_block(tprfile, xtcfile, blockslice, func, *args, **kwargs):
    """
    Load a trajectory block and apply an analysis function.
    Returns the result for that block.

    Arguments:
        tprfile    : Path to GROMACS .tpr file.
        xtcfile    : Path to trajectory (.xtc) file.
        blockslice : Range/slice of frames to analyze.
        func       : Analysis function to apply.
        *args      : Additional positional arguments for func.
        **kwargs   : Additional keyword arguments for func.

    Returns:
        result     : Output of func for the given trajectory block.
    """
    u_block = Universe(tprfile, xtcfile )
    block_start = max(0, blockslice.start)
    block_stop = min(u_block.trajectory.n_frames, blockslice.stop)
    print(f"Analyzing block: {block_start} to {block_stop}")
    return func(u_block, block_start, block_stop, *args, **kwargs)

def split_traj(start, stop, n_jobs):
    """
    Split trajectory frame range into blocks for parallel processing.
    Divides frames evenly across n_jobs workers.
    Returns a list of frame ranges.

    Arguments:
        start   : Starting frame index.
        stop    : Ending frame index.
        n_jobs  : Number of blocks/workers.

    Returns:
        blocks  : List of range objects defining frame blocks.
    """
    n_frames = stop - start
    n_blocks = n_jobs 
    n_frames_per_block = n_frames // n_blocks
    blocks = [range(i * n_frames_per_block, (i + 1) * n_frames_per_block) for i in range(n_blocks-1)]
    blocks.append(range((n_blocks - 1) * n_frames_per_block, n_frames))
    return blocks

def an_perblock(logger, tprfile, xtcfile, n_jobs, func, start=0, stop=None, *args, **kwargs):
    """
    Run trajectory analysis in parallel using Dask.
    Splits trajectory into blocks and applies a function per block.
    Combines results into a single array.

    Arguments:
        logger  : Logger object for debug output.
        tprfile : Path to GROMACS .tpr file.
        xtcfile : Path to trajectory (.xtc) file.
        n_jobs  : Number of parallel workers.
        func    : Analysis function to apply per block.
        start   : Starting frame index (default: 0).
        stop    : Ending frame index (default: full trajectory).
        *args   : Additional positional arguments for func.
        **kwargs: Additional keyword arguments for func.

    Returns:
        results : NumPy array containing concatenated results from all blocks.
    """
    logger.debug(f"Calling an_perblock. Passing files: TPR: {tprfile}, XTC: {xtcfile}")
    setup_dask(n_jobs)
    logger.debug(f"Dask setup complete.")
    u_temp = Universe(tprfile, xtcfile) 
    logger.debug(f"u_temp Universe loaded.")
    if start is None or start < 0:
        start = 0
        logger.debug(f"Start frame not provided or negative. Using start frame: {start}")
    if stop is None: 
        stop = u_temp.trajectory.n_frames
        logger.debug(f"Stop frame not provided. Using total number of frames: {stop}")
    blocks = split_traj(int(start), int(stop), int(n_jobs)) 
    logger.debug(f"Trajectory split into {len(blocks)} blocks for parallel processing.")
    logger.debug(f"Blocks: {blocks}")
    if n_jobs == 1:
        logger.debug("Running single-threaded for debugging")
        result = analyze_block(tprfile, xtcfile, blocks[0], func, *args, **kwargs).compute()
        return np.array([result])
    jobs = [analyze_block(tprfile, xtcfile, bs, func, *args, **kwargs) for bs in blocks]
    logger.debug(f"Dask jobs created for each block.")
    results = dask.compute(*jobs)
    logger.debug(f"Dask computation complete. Number of results: {len(results)}")
    return np.concatenate(results)
