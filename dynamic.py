import sys
import time

for i in range(10):
    sys.stdout.write("\r{0}>".format("="*i))
    print(i)
    sys.stdout.flush()
