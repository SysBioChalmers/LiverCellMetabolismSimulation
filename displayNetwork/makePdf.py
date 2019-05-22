import os
import sys
from graphviz import Source
os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'

path = sys.argv[1]
fileEnding = sys.argv[2]
file = open(path, 'r')#READING DOT FILE
text=file.read()
s = Source(text, filename=path, format=fileEnding)

s.view()
file.close()
