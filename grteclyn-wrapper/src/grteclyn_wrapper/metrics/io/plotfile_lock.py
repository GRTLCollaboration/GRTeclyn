"""Thread lock for yt plotfile reads during parallel GPU evaluation."""

from __future__ import annotations

import threading

# yt is not thread-safe; the optimizer evaluates a generation across GPUs using
# threads, so concurrent yt.load() calls in different threads race and the
# evolved-FTL/effective-EC reads silently fail.  Serialize just the plotfile
# (yt) reads with this lock so the GPU simulations still run fully in parallel.
PLOTFILE_READ_LOCK = threading.Lock()
