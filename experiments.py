# import python packages
import pickle
from deap import base, creator, gp, tools
from deap import algorithms as algo
from random import random, seed, shuffle
import random as rand
import networkx as nx
import numpy as numpy
import copy as copy
import operator
import math as math
from sets import Set
from joblib import Parallel, delayed
import argparse as argparse
import scipy.stats as stat
from sets import Set
import matplotlib.pyplot as plt
import requests

import sklearn.mixture as mix
# import other pieces of our software
import networkConstructor as nc
import utils as utils
import simulation as sim
import murphyReader as mr
import motif_cutter as mc
import pathway_analysis as pa
import GA as ga
def deltaTester(ultimate, i, model,  newSSS, params, minny):
	print(i)
	modifyFlag=False
	copied=list(ultimate)
	copied[i]=1-copied[i]
	for node in range(0,len(model.nodeList)):
		better=True
		#get start and end indices for node in individual
		if node==(len(model.nodeList)-1):
			end=len(ultimate)-1
		else:
			end=model.individualParse[node+1]
		start=model.individualParse[node]
		if model.andLenList[node]>0:
			if sum(copied[start:end])>0:
				newtot=numpy.sum(evaluateByNode(copied, params.cells, model,  newSSS, params, KOlist, KIlist))
				if newtot<minny:
					modifyFlag=True
	return modifyFlag

def simTester(graph, name):
	#creates a model, runs simulations, then tests reverse engineering capabilities of models in a single function
	#trials is the number of different individuals to try
	#samples is the number of different initial conditions to provide per trial
	#graph specifies the network we are testing. 

	# set up params, sss, model, samples


	params=sim.paramClass()
	# rewire graph if necessary
	if params.rewire:
		graph=nc.rewireNetwork(graph)

	sss=utils.synthesizeInputs(graph,params.samples)
	model=sim.modelClass(graph,sss)
	samples=params.samples

	# set up necessary lists for output
	truthList= [] 
	testList =[]
	devList=[]

	knockoutLists=[]
	knockinLists=[]

	for q in range(samples):
		temRand=rand.randint(0,len(model.nodeList))
		knockoutLists.append([])
		knockinLists.append([])


	# loop over number of times we want to generate fake data and perform sequence of events
	for i in range(0,params.trials):
		#generate random set of logic rules to start with
		individual=genBits(model)
		for node in range(0,len(model.nodeList)):
			if model.andLenList[node]>0:
				if node==len(model.nodeList)-1:
					end=len(model.nodeList)
				else:
					end=model.individualParse[node+1]
				if numpy.sum(individual[model.individualParse[node]:end])==0:
					individual[model.individualParse[node]]=1

		# generate Boolean model for this trial
		output=runProbabilityBooleanSims(individual, model, samples, params.cells, params, knockoutLists, knockinLists)
		# copy new output into newSSS and initial values
		newSSS=[]
		for k in range(0,samples):
			newSS=copy.deepcopy(sss[k])
			for j in range(0,len(model.nodeList)):
				newSS[model.nodeList[j]]=output[k][j]
			newSSS.append(newSS)
		newInitValueList=[]
		for j in range(0,len(sss)):
			newInitValueList.append([])
		for j in range(0,len(model.nodeList)):
			for k in range(0,len(sss)):
				ss=newSSS[k]
				if  model.nodeList[j] in sss[0]:
					newInitValueList[k].append(ss[model.nodeList[j]])
				else:
					newInitValueList[k].append(0.5)
		model.initValueList=newInitValueList
		#perform GA run
		dev, bruteOut=GAsearchModel(model, newSSS, params, knockoutLists, knockinLists)
		# add output of this trial to lists
		truthList.append(individual)
		testList.append(list(bruteOut))
		devList.append(dev)
	# set up output and save as a pickle
	outputList=[truthList,testList, devList,model.size, model.nodeList, model.individualParse,model.andNodeList ,model.andNodeInvertList, model.andLenList,	model.nodeList, model.nodeDict, model.initValueList]
	pickle.dump( outputList, open( name+"_output.pickle", "wb" ) )

def partitionTester(graph,name,ssDict):
	# first trim the graph to just what is in the dataset and the dataset to just what is in the graph

	if len(graph.nodes())>1:	
		dev3, bruteOut3, model3= partitionTest(graph, ssDict) # divide by max
		pickle.dump( [[dev3],[bruteOut3],[model3]], open( name+"_output.pickle", "wb" ) )
	else:
		print('not enough overlap')

def partitionTest(graph,  ssDict):

	# ssDict, valueLister, ones, zeros= mr.constructOmicsInput(geneDict, maxRem, divider)
	params=sim.paramClass()
	
	graph=nc.simplifyNetwork(graph, ssDict) # simplify the graph
	graph=mc.cutMotifs(graph)
	print(graph.nodes())
	keyList=ssDict.keys()
	# generate SSS which is a list of SSs where SS is a dict pointing to some values
	sss=[{} for i in range(len(ssDict[keyList[0]]))]
	newInitValueList=[[] for i in range(len(ssDict[keyList[0]]))]
	for i in range(len(sss)):
		for key in keyList:
			if key in graph.nodes():
				sss[i][key]=ssDict[key][i]
	
	model=sim.modelClass(graph,sss)
	# generate empty knockout and knockin lists for support
	knockoutLists=[]
	knockinLists=[]
	for q in range(len(sss)):
		knockoutLists.append([])
		knockinLists.append([])

	# generate a set of initial values which are just the values from the SS in order. 
	newInitValueList=[]
	for j in range(0,len(sss)):
		newInitValueList.append([])
	for j in range(0,len(model.nodeList)):
		for k in range(0,len(sss)):
			ss=sss[k]
			if  model.nodeList[j] in sss[0]:
				newInitValueList[k].append(ss[model.nodeList[j]])
			else:
				newInitValueList[k].append(0.5)
	model.initValueList=newInitValueList
	print('setup successful')
	#perform GA run
	dev, bruteOut= ga.GAsearchModel(model, sss, params, knockoutLists, knockinLists)
	return dev, bruteOut, model

def findPathways(geneDict):
	returnerlist=[]
	aliasDict={}
	dict1={}
	nc.parseKEGGdicthsa('inputData/hsa00001.keg',aliasDict,dict1)
	dict2={}
	nc.parseKEGGdict('inputData/ko00001.keg',aliasDict,dict2)
	namelist=pa.find_overlaps('filtered.c2.cp.kegg.v3.0.symbols.gmt',geneDict)
	print('num of overlap nodes')
	print(len(namelist))
	for name in namelist:
		returnerlist.append(pa.retrieveGraph(name,aliasDict,dict1,dict2, geneDict))
	return returnerlist

def constructBinInput(filename):
	binput = open(filename, 'r')
        ssDict = {}
	lines = binput.readlines()
	for line in lines:
		geneCount = line.split('\t')
		#i = len(geneCount) - 1
		#ssDict[geneCount[0]] = geneCount[1:i]
		ssDict[geneCount[0]] = geneCount[1:-1]
	return(ssDict)

if __name__ == '__main__':
	import time
	start_time = time.time()

	# use the command line to input the file that you want
	#Is this right?
	#python experiments.py graphname(to be given) binData\kMeans.bin iterNum(to be given)

	parser = argparse.ArgumentParser()
	parser.add_argument("graph")
	parser.add_argument("datafile")
	parser.add_argument("iterNum")

	results = parser.parse_args()
	graphName=results.graph
	fileName=results.datafile
	iterNum=int(results.iterNum)
	
	name=graphName[:-8]+fileName+'_'+results.iterNum ##

	# insert your code to get a ssDict here.
	
	# use fileName to read in the discretized file... you want to standardize the format
	#ssDict= SOMETHING
	ssDict = constructBinInput(fileName)
	
	# edit the below name to include the pertinent piece of the filename of discretized data
	#outname=graphName[:-8]+'_'+results.iterNum
	outname=graphName[:-8]+fileName+'_'+results.iterNum
	print(name) #
	graph = nx.read_gpickle(name)
	print(len(graph.nodes()))
	partitionTester(graph,outname,ssDict)
	print("--- %s seconds ---" % (time.time() - start_time))
