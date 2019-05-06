import numpy as np
import sys
x = np.loadtxt(sys.argv[1],delimiter=',', dtype=int)
y = np.dot(x,x.T)
np.savetxt(sys.argv[2],y, fmt='%d')