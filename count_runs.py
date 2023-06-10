#!/usr/bin/env python3

import sys, time, argparse, subprocess, os, datetime

Description = """
Tool to count n/r of a bwt.
"""

def noRuns(filepath):

    if not os.path.isfile(filepath):
        print("File path {} does not exist. Exiting...".format(filepath))
        sys.exit()
  
    r = 0
    with open(filepath, 'rb') as fp:
        r = 1
        char = fp.read(1)           
        while 1:
            prev = char
            char = fp.read(1)           
            
            if not char:  
                break

            if prev != char:
                r +=1
                prev = char

    return r
