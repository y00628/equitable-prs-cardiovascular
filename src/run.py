import pandas as pd
import numpy as np
import os
import subprocess
import urllib.request
import data
import compute_weights_win
import extract_data
import test
import sys

if __name__ == '__main__':
    args = sys.argv[1:]

    if 'test' in args:
        test.test()
