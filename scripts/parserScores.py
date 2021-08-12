#!/usr/bin/env python2

from __future__ import print_function
import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import sys
import csv

from optparse import OptionParser


class CmdLineParser(object):
	def __init__(self):
		self.parser = OptionParser(usage='usage: python2 parserScores.py scores.csv outputfile.csv')

	def addOption(self, *args, **kwargs):
		self.parser.add_option(*args, **kwargs)

	def parseArgs(self):
		(options, args) = self.parser.parse_args()
		options.input = args[0]
		options.output = args[1]
		return options


def parseScores(fp, outputfile):
	"""
	read a CSV of scores
	the format is as follows:
	instance, seed, model, score, time, run_number
	"""
	firstRow = [0] # n_models, model1, model2, ...
	instances = [] # name of the instances(with the seed)
	runs = []	# rows 


	for row in fp:
		row = row.strip().split(', ')
		instance = row[0] + "_" + row[5]
		model = row[2]
		time = row[4]

		try: #search for the model's index
			model_index = firstRow.index(model)-1
		except ValueError:
			#print("new model found", model)
			firstRow.append(model)
			firstRow[0] += 1
			model_index = firstRow[0]-1

		try: #search for the instance's index
			instance_index = instances.index(instance)
			l = len(runs[instance_index])
			if (model_index >= l): 
				# add extra cells
				runs[instance_index].extend([None]* (model_index-l+1))
			runs[instance_index][model_index] = time
		except ValueError:
			#print("new instance:", instance)
			instances.append(instance)
			instance_index = instances.index(instance)	
			runs.append([time])

	# create table
	with open(outputfile, mode='a') as file:
		writer = csv.writer(file, delimiter=',')
		writer.writerow(firstRow)
		# write runs
		for i in range(len(instances)):
			writer.writerow([instances[i]] + runs[i])
		

def main():
	parser = CmdLineParser()
	opt = parser.parseArgs()
	parseScores(open(opt.input, 'r'), opt.output)


if __name__ == '__main__':
	main()
