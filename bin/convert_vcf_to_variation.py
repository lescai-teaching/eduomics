#!/bin/env python

import sys
import random
import argparse

def get_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", help="vcf file to transform variants",
						action="store", dest = "input")
	parser.add_argument("-o", "--output", help="output reformatted variation to file - can be existing one",
						action="store", dest = "output")
	parser.add_argument("-n", "--normal", help="choose if writing lines for normal individuals or not",
						action="store", dest = "normal")
	args = parser.parse_args()
	return(args)

def parse_variant_info(string):

	fields = string.split("\t")
	chr = fields[0]
	pos = fields[1]
	ref = fields[3]
	alt = fields[4]
	genos = ["het","het","het","homo"]
	genonormal = random.choice(genos)
	genodisease = random.choice(genos)
	descriptor = ""
	type = ""
	normline = ""
	disline = ""

	if len(ref) > len(alt):
		type = "d"
		pos = str(int(pos) + 1)
		descriptor = str( len(ref) - len(alt) )
		normline = "\t".join([ type, "normal", chr, pos, descriptor, genonormal])
		disline = "\t".join([ type, "disease", chr, pos, descriptor, genodisease])
	elif len(ref) < len(alt):
		type = "i"
		descriptor = alt.removeprefix(ref)
		normline = "\t".join([ type, "normal", chr, pos, descriptor, genonormal])
		disline = "\t".join([ type, "disease", chr, pos, descriptor, genodisease])
	elif len(ref) == len(alt):
		type = "s"
		normline = "\t".join([ type, "normal", chr, pos, ref, alt, genonormal])
		disline = "\t".join([ type, "disease", chr, pos, ref, alt, genodisease])
	return(normline, disline)



def main():
	args = get_arguments()
	filein = open(args.input)
	fileout = open(args.output, "a")
	for line in filein:
		line = line.strip()
		if line.startswith("#"):
			continue  # Skip header lines
		normline, disline = parse_variant_info(line)
		if "." in disline:
			continue
		if args.normal:
			fileout.write(normline + "\n")
			fileout.write(disline + "\n")
		else:
			fileout.write(disline + "\n")
	filein.close()
	fileout.close()


if __name__ == '__main__':
	main()

