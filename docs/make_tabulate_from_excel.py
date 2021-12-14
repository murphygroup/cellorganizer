from tabulate import tabulate
import pandas 
import sys
from os import getcwd

inpt = sys.argv[1]

#rapido y sucio
try:
	sheet = sys.argv[2]
except:
	sheet = 'Sheet1'

try:
	fill_na_string = sys.argv[3]
except:
	fill_na_string = 'N/A'

df = pandas.read_excel(inpt, sheet)
df = df.fillna( fill_na_string )
df = df.replace( False, '' )
df = df.replace( True, 'True' )
body = df.values.tolist()
header = df.columns.values.tolist()

print(tabulate(body, header, tablefmt="grid"))
