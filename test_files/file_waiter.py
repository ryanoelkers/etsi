# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 00:21:25 2022

@author: Jimmy Ardoin
"""

import os
import time
from config import Configuration
from scripts.clean import Clean


def file_wait(num_files=1, path=Configuration.RAW_DIRECTORY, timeout=0):
    """
    Will wait for timeout minutes for num_files to appear in path. Then cleans
    those files with the functions from clean.py. If timeout is 0 or negative
    will wait indefinitely.

    Parameters
    ----------
    num_files : int, optional
        number of files to wait on. The default is 1
    path : str, optional
        Path of the directory to wait in. The default is the raw directory
        found in config.py.
    timeout : float, optional
        How many minutes to wait before timeing out. Any number <= 0 will make
        it run indefinitely. The default is 0.

    Returns
    -------
    None. Will break when timeout minutes have been reached (if specified) or
    when there are num_files in path.

    """

    start_time = time.time()

    while True:
        # Get contents of the directory and removes this script from the list
        # if no path is specified
        dr = os.listdir(path)
        if path == '.':
            dr.remove(os.path.basename(__file__))

        # If there are num_files there, it will print the contents and break
        if len(dr) >= num_files:
            Clean.clean_images(dark_subtract='Y', sky_subtract='Y')
            return

        # If there is no file then it will wait a little and then resume the
        # loop
        else:
            time.sleep(1)

        # Check current time to see if timeout minutes has passed
        if timeout > 0:
            ctime = time.time()
            if ctime > (start_time + timeout * 60):
                return

