import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def main():
    pathRes1 = sys.argv[1]
    pathRes2 = sys.argv[2]

    try:
        res1 = np.loadtxt(pathRes1)
        res2 = np.loadtxt(pathRes2)

        assert len(res1) == len(res2)
        yAbs = [xi-xj for xi, xj in zip(res1, res2)]

        absoluteFault = max(map(abs,yAbs))
        print("%.15f" % absoluteFault)

    except AssertionError:
        print ('ERROR! Не совпадают размерности! ' + str(len(res1)) + "!=" + str(len(res2)))

if __name__ == '__main__':
    main()
