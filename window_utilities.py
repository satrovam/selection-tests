#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Zde mam ulezeno par funkci z ruznych scikit-allel sourse codes,
ktere se obcas hodi.

"""

# This is the function allel uses to compute windows
# I will use it to compute the respective window positions
def index_windows(values, size, start, stop, step):
    """Convenience function to construct windows for the
       :func:`moving_statistic` function.

    """

    # determine step
    if stop is None:
        stop = len(values)
    if step is None:
        # non-overlapping
        step = size

    # iterate over windows
    for window_start in range(start, stop, step):

        window_stop = window_start + size
        if window_stop > stop:
            # ensure all windows are equal sized
            return

        yield (window_start, window_stop)