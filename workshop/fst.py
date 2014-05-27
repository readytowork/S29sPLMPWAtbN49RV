import os
os.system("date")

f = os.popen('date')
now = f.read()
print "Today is ", now	

import subprocess
subprocess.call("ls")
subprocess.call(["ls", "-alh", "."])

subprocess.call(["./run.sh", "exp08.cpp"])

p = subprocess.Popen(["ls", "-alh", ".."], stdout=subprocess.PIPE)
output, err = p.communicate()
print "*** Running ls -l command ***\n", output
