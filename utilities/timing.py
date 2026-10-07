"""Lightweight wall-clock profiling for pipeline steps.
"""

import time
from contextlib import contextmanager

from utilities import log


@contextmanager
def timed(label: str):
    t0 = time.perf_counter()
    try:
        yield
    finally:
        log.info(f"  [timing] {label}: {time.perf_counter() - t0:.2f}s")
