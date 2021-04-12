import random
import os
import sys


def generate():
    directory = "../data/random"

    # Check whether the directory exists
    if not os.path.exists(directory):
        os.makedirs(directory)

    n = int(sys.argv[1])

    size = [20, 30, 40, 50, 60, 70, 80, 90, 100, 110]

    a = random.randint(-10000, -1)
    b = random.randint(1, 10000)

    for s in size:

        subdirectory = directory + "/" + str(s)

        # Check whether the subdirectory exists
        if not os.path.exists(subdirectory):
            os.makedirs(subdirectory)

        for i in range(n):

            filename = "instance" + str(i + 1)
            path = subdirectory + "/" + filename + ".tsp"

            if not os.path.isfile(path):
                file = open(path, "x")
            else:
                file = open(path, "w")

            file.write("NAME : " + filename + "\n")
            file.write("COMMENT : " + str(s) + " random nodes\n")
            file.write("TYPE : TSP\n")
            file.write("DIMENSION : " + str(s) + "\n")
            file.write("EDGE_WEIGHT_TYPE : EUC_2D\n")
            file.write("NODE_COORD_SECTION\n")

            for k in range(0, s):
                x = random.randint(a, b)
                y = random.randint(a, b)
                file.write(str(k + 1) + " " + str(x) + " " + str(y) + "\n")

            file.write("EOF")
            file.close()


if __name__ == '__main__':
    generate()
