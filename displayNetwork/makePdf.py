import os
from graphviz import Source
os.chdir('D:\stories\Glutamine\matlab\displayNetwork')
file = open('test.dot', 'r')#READING DOT FILE
text=file.read()
s = Source(text, filename='test.dot', format="pdf")

s.view()
file.close()
