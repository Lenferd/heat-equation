# Скрипт удаляет промежуточные значения (идет с шагом 2, а не 1),
# таким образом позволяя сгенерировать файл с большей точностью для валидации
# файлов с меньшей точностью (большая точность шага по пространству)

# Задаваемые параметры:
# nx текущее значение количества точек в входном файле

import numpy as np
import sys

if __name__ == '__main__':
    # Set this parametrs
    nx = 17
    infile = "function2-129m17.txt"
    outfile = "function2-129m9.txt"

    if len(sys.argv) != 4:
        print("Using default numbers. Please take a look inside script")
        print("Or use $ python3 delscript.py <nx> <infile name> <outfile name>\n")
    else:
        nx = int(sys.argv[1])
        infile = sys.argv[2]
        outfile = sys.argv[3]

    print("Nx: " + str(nx) +" infile: " + infile + " outfile: " + outfile)
    data = np.fromfile(infile, float, sep='\n',)
    data.shape = (nx, nx, nx)

    result_data = data[::2, ::2, ::2]

    np.ndarray.tofile(result_data, outfile, sep='\n', format="%.15e")
