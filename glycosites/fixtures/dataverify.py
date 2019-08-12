import re

filename = "HumanAll.txt"
f = open(filename,'r')
for line in f:
    i = re.split('\t',line.rstrip())
    x = re.split(';',i[8])
    y = re.split(';',i[9])
    if len(x)!= len(y):
        print(i[1],len(x),len(y),i[8],i[9])