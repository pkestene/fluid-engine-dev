#!/usr/bin/env bash

mkdir build
cd build
cmake .. -DENABLE_CUDA=OFF
make -j`nproc`
bin/unit_tests
cd ..
pip install --user -r requirements.txt
pip install --user .
python src/tests/python_tests/main.py
