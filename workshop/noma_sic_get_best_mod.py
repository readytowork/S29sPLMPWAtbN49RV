#!/usr/bin/python
import os
import subprocess

f = os.popen('date')
now = f.read()
print "Today is ", now	

subprocess.call(["./run.sh", "exp09.cpp"])
n_bits = 32768
for dist_1 in range(400, 2000, 100):
	for dist_2 in range(dist_1, 2000, 100):
		for mod_1 in range(1, 3, 1):
			for mod_2 in range(1, 3, 1):
				for rate_1 in range(1, 3, 1):
					for rate_2 in range(1, 3, 1):
						subprocess.call(["./exp09", str(n_bits), str(dist_1), str(dist_2), str(mod_1), str(mod_2), str(rate_1), str(rate_2)])
