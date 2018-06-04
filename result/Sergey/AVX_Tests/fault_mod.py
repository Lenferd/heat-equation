import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def main():
    if len(sys.argv) != 4:
        print("Не верное количество параметров!")
        print("settings.txt file1.txt file2.txt")
        return

    settingPath = sys.argv[1]
    with open(settingPath, 'r') as file:
        pattern = re.compile('[A-Za-z]+=-?\d+')
        setting = { line.split('=')[0] : float(line.split('=')[1])
                    for line in file if pattern.match(line) }

    xStart = setting['XSTART']
    xFinish = setting['XEND']
    NX = setting['NX']

    x = np.linspace(xStart, xFinish, NX)

    pathRes1 = sys.argv[2]
    pathRes2 = sys.argv[3]

    try:
        res1 = np.loadtxt(pathRes1)
        res2 = np.loadtxt(pathRes2)

        assert len(res1) == len(res2)

        yAbs = [xi-xj for xi, xj in zip(res1, res2)]
        absoluteFault = max(map(abs,yAbs))

        yRelat = [(xi - xj) / max(xi, xj) for xi, xj in zip(res1, res2)]
        # yRelatAll = [(xi, xj, (xi-xj)/max(xi, xj)) for xi, xj in zip(res1, res2)]
        relativeFault = max(map(abs,yRelat))
        # relativeFault = abs, lambdaL
        # relativeFaultMax = max(yRelat)
        # relativeFaultMin = min(yRelat)

        # print(relativeFault)
        value = [(xi, xj) for xi, xj in zip(res1, res2) if abs((xi - xj) / max(xi, xj)) == relativeFault]
        # value = [ for xi, xj in zip(res1, res2)

        print("абсолютная:\t%.15f" % absoluteFault)
        print("относительная:\t%.15f При\t%s, %s" % (relativeFault , value[0][0] , value[0][1]))


        ####################################################################
        #                    Рисование графиков                            #
        ####################################################################

        assert len(x) == len(yAbs)

        plt.figure(num = 'FAULT', facecolor = (1, 1, .54))

        plt.subplot(221)
        plt.plot(x, yAbs, label ='absolut', color = 'green')

        plt.legend(loc = 2)
        plt.xlabel('x', fontsize = 14)
        plt.grid(True)
        ax = plt.gca()
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{:.0e}'.format(x)))
        # plt.yscale('log')

        plt.subplot(222)
        plt.plot(x, yRelat, label = 'relative', color = 'red')
        plt.legend(loc = 2)
        plt.xlabel('x', fontsize=14)
        plt.grid(True)

        plt.subplot(223)
        plt.plot(x, res1, label = sys.argv[1], color = 'blue')
        plt.legend(loc = 2)
        plt.xlabel('x', fontsize = 14)
        plt.grid(True)

        plt.subplot(224)
        plt.plot(x, res2, label = sys.argv[2], color = 'brown')
        plt.legend(loc = 2)
        plt.xlabel('x', fontsize = 14)
        plt.grid(True)

        plt.show()
    except AssertionError:
        print ('ERROR! Не совпадают размерности! ' + str(len(res1)) + "!=" + str(len(res2)))


if __name__ == '__main__':
    main()
