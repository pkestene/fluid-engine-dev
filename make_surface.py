from __future__ import print_function
import os

EXE_CMD = "bin/particles2obj"
POS_DIR = "/Users/doyubkim/Desktop/example6_res128/"
OBJ_DIR = "/Users/doyubkim/Desktop/example6_res128_obj/"
RESX = 300
LENX = 1.0
GRID_SIZE = LENX / RESX
METHODS = ["anisotropic", "sph", "spherical", "zhu_bridson"]
KERNELS = [0.015, 0.015, 0.0075, 0.025]
OPTIONS = ["anisotropic,0.6,0.9", "sph,0.5", "spherical", "zhu_bridson,0.23"]


def exp(m):
    for i in range(0, 240):
        posname = "%sframe_%06d.pos" % (POS_DIR, i)
        objname = "%sframe_%06d.%s.obj" % (OBJ_DIR, i, METHODS[m])

        cmd = "%s -i %s -o %s -r %d -g %f -k %f -m %s" % (
            EXE_CMD, posname, objname, RESX, GRID_SIZE, KERNELS[m], OPTIONS[m])
        os.system(cmd)


def main():
    if not os.path.exists(OBJ_DIR):
        os.makedirs(OBJ_DIR)

    for i in range(1, 4):
        exp(i)


if __name__ == "__main__":
    main()
