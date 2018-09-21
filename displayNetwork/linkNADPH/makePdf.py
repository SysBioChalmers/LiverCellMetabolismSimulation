import os
from graphviz import Source
os.chdir('D:\stories\Glutamine\matlab\displayNetwork\linkNADPH')
file = open('test.dot', 'r')#READING DOT FILE
text=file.read()
s = Source(text, filename='test.dot', format="pdf")
s.view()
file.close()


file = open('test_legend.dot', 'r')#READING DOT FILE
text=file.read()
s = Source(text, filename='test_legend.dot', format="pdf")
s.view()
file.close()