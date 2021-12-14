#!/usr/bin/python

from tabulate import tabulate
import sys

if sys.argv[1] == '--help':
	print 'Usage: python make_tabulate.py [input]'
	print 'Please add proper HEADER to your text file if you need them'
	print 'Please use SINGLE SPACE as separator'
	sys.exit()

inpt = sys.argv[1]
inpt = open(inpt).read().splitlines()

header = inpt.pop(0).split(' ')
table = []
for i in inpt:
	i = i.split(' ')
	table.append(i)

print tabulate(table, header, tablefmt="grid")