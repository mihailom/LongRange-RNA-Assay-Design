#for generating list of unpaired base pairing probabilities from NUPACK 
from collections import OrderedDict
import argparse

parser = argparse.ArgumentParser(description='clean nupack .ppairs output')
parser.add_argument('-i',
	action='store',
	dest = 'i',
	required=True,
	type=str,
	help='nupack .ppairs output file')

#parser.add_argument('-o',
#	action='store',
#	required=True,
#	type=str,
#	help='output file of cleaned nupack .clean.txt')

options = parser.parse_args()

file = open(options.i,'r') #Set file to match RNA query

linecounter = 0

for line in file:
    linecounter += 1
    newline = line.strip('\n')
    # print(newline)
    newline_list = newline.split('\t')
    if linecounter == 15:
        RNA_length = int(newline_list[0])
        unpair_prob_dict = OrderedDict()
        # print (RNA_length)
    if linecounter >= 16 and newline_list[1] == str(RNA_length+1):
        unpair_prob_dict[newline_list[0]] = float(newline_list[2])

file.close()
#print(unpair_prob_dict)
#print ('\n')

root_outfile_name = options.i[0:len(options.i)-7]+'.clean.txt'

outfile = open(root_outfile_name, "w")

for i in range(1,RNA_length + 1):
#    print i
    if unpair_prob_dict.has_key(str(i)):
        keystring = str(i)
#        print ('Yes! key is present')
        outfile.write(keystring + '\t' + str(round(unpair_prob_dict[keystring], 3)) + '\n')
    else:
        outfile.write(str(i) + '\t' + '0.0' + '\n')
        
outfile.close()
    
                  
