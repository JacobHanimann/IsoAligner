import sys
print(sys.path)
sys.path.insert(0,'/Users/jacob/PycharmProjects/IsoAligner/IsoAligner_core')
print(sys.path)
import os

# get current directory
path = os.getcwd()
print("Current Directory", path)

# prints parent directory
print(os.path.abspath(os.path.join(path, os.pardir)))
print(os.pardir)
print(os.path.abspath("hello"))
print(os.path.join(path,os.pardir))
print(os.path.dirname(path))

print(os.curdir)

print(os.path.abspath(os.curdir))
os.chdir("../../../IsoAligner_core")
print(os.path.abspath(os.curdir))