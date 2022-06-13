import os
import sys


print("--- SOME DIAGNOSTIC VALUES: ---")

print("\nsys.executable: ", sys.executable)

print("\ncwd (Current Working Directory): ", os.getcwd())

print("\nsys.path: ", sys.path)
print("\nsys.path again, this time on separate lines: ")
for p in sys.path:
    print("    ", p)
