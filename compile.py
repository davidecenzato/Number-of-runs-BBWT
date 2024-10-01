#!/usr/bin/env python3

import sys, subprocess, os

os.chdir("cais")

os.chdir("external/sdsl-lite")
subprocess.call("mkdir installed".split())
subprocess.call("./install.sh installed/".split())

os.chdir("../../")
subprocess.call("cp ../internal/main.cpp ./".split())
subprocess.call("cp ../internal/Makefile ./".split())
subprocess.call("make".split())

os.chdir("../")
subprocess.call("mkdir outdir".split())