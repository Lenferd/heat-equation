import numpy as np
from math import pi, sin
import os


def function(x, y, z, xmax, ymax, zmax):
    return sin(2 * pi * (x - xmax) / xmax) + sin(2 * pi * (y - ymax) / ymax) + sin(2 * pi * (z - zmax) / zmax)


def main():
    #pathSetting = os.path.join(os.pardir, "initial_for_tests", "setting5.ini")
    pathSetting = "setting2.ini"

    with open(pathSetting) as file:
        setting = {line.split('=')[0]: float(line.split('=')[1]) for line in file}

    xStart = setting['XSTART']
    xEnd = setting['XEND']
    NX = setting['NX']

    yStart = setting['YSTART']
    yEnd = setting['YEND']
    NY = setting['NY']

    zStart = setting['ZSTART']
    zEnd = setting['ZEND']
    NZ = setting['NZ']

    X = np.linspace(xStart, xEnd, NX, dtype=float)
    Y = np.linspace(yStart, yEnd, NY, dtype=float)
    Z = np.linspace(zStart, zEnd, NZ, dtype=float)

    U = np.array([function(x, y, z, xEnd, yEnd, zEnd) for z in Z for y in Y for x in X])

    #pathOutput = os.path.join(os.pardir, "initial_for_tests", 'function5.txt')
    pathOutput = 'function2.txt'
    np.savetxt(pathOutput, U, fmt='%.15e', delimiter='\n')


if __name__ == "__main__":
    print("Start")
    main()
print("Finish")
