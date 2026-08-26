"""quant209: PCB GC-MS/MS quantification pipeline.

Pin BLAS/OpenMP thread counts to 1 before anything imports numpy. The host
is memory-constrained and multi-threaded BLAS causes contention/OOM. This
module must stay import-cheap (no other imports) so every entry point can
`import quant209` first without cost.
"""

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

__version__ = "1.0.0"
