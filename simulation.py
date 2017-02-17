#import from other parts of ruleMaker
import utils as utils
#import python modules
from random import random
from random import shuffle
import operator
import networkx as nx
import itertools as itertool

class modelClass:
     def __init__(self,graph, sss): 

     	#remove self loops from the graph
		for node in graph.nodes():
			repeat=True
			while(repeat):
				repeat=False
				if node in graph.successors(node):
					graph.remove_edge(node,node)
					repeat=True
		
		#set up empty lists and dicts for later
		evaluateNodes=[] #list of nodes that need to be compared to the steady state values (sss)
		individualParse=[] # list of the number of shadow and nodes that contribute to each node, in order by index num
		earlyEvalNodes=[] # nodes that don't have initial value and need to be re-evaluated early on in the simulation.. 
		andNodeList=[] #a list of the shadow nodes that represent and relations between incoming edge
		andNodeInvertList=[] # keeps track of which incoming nodes for each node need to be inverted
		andLenList=[] # keeps track of how many nodes are coming into each shadow AND node
		nodeList=graph.nodes()#define the node list simply as the nodes in the graph. 
		nodeDict={} #identifies names of nodes with their index in the node list- provide name, get index
		for i in range(0,len(nodeList)):
			nodeDict[nodeList[i]]=i #constructs the node dict so we can easily look up nodes
		counter=int(0) #keeps track of where we are in the generic individual
		initValueList=[] #starting states for nodes
		for j in range(0,len(sss)): #construct empty lists to stick values in later for intiial value list
			initValueList.append([])
		
		#find all possible combinations of upstream contributors for each node. These become the shadow And nodes
		for i in range(0,len(nodeList)):
			preds=graph.predecessors(nodeList[i]) # get predecessors of node. 
			if len(preds)>15: #handle case where there are too many predecessors by truncation
				preds=preds[1:15]
			for j in range(0,len(preds)):
				preds[j]=nodeDict[preds[j]]
			# the followign section constructs a list of possible node orders
			# this is accomblished by finding all possible subsets of the list of predecessor nodes
			withNones = zip(preds, itertool.repeat('empty'))
			possibilities=list(itertool.product(*withNones))
			for j in range(0,len(possibilities)):
				possibilities[j]=list(possibilities[j])
				while 'empty' in possibilities[j]:
					possibilities[j].remove('empty')
				while [] in possibilities[j]:
					possibilities[j].remove([])
			while [] in possibilities:
				possibilities.remove([])
			
			# create a list of the activities of each node and store alongside the contributors to each and node for easy reference later
			activities=[] #list to store activities of nodes (a vs i)
			for sequence in possibilities:
				activity=[]
				for node in sequence:
					if graph.edge[nodeList[node]][nodeList[i]]['signal']=='a':
						activity.append(False)
					else:
						activity.append(True)
				activities.append(activity)
			andNodeList.append(possibilities)
			andNodeInvertList.append(activities)
			andLenList.append(len(possibilities))
			
			# construct the list of lengths of possibilties for each node, add to the counter that keeps track of how many bits are necessary
			individualParse.append(counter)
			counter=counter+len(possibilities)
			if  nodeList[i] in sss[0].keys():
				evaluateNodes.append(i)
			else:
				earlyEvalNodes.append(i)
		self.size=counter
		self.evaluateNodes=evaluateNodes #list of nodes that need to be compared to the steady state values (sss)
		self.individualParse=individualParse #index of start value of current node on the individual
		self.andNodeList=andNodeList # shadow add node inputs
		self.andNodeInvertList=andNodeInvertList # keeps track of which incoming nodes for each node need to be inverted
		self.andLenList=andLenList # keeps track of length of above inputOrderList for each node
		self.earlyEvalNodes=earlyEvalNodes
		self.nodeList=nodeList #define the node list simply as the nodes in the graph. 
		self.nodeDict=nodeDict #identifies names of nodes with their index in the node list.. provide name, get index
		self.initValueList=initValueList #puts an empty and correctly structured initValueList together for later population. 

class paramClass:
	def __init__(self):     
		self.simSteps=100 # number of steps each individual is run when evaluating
		self.generations=15 # generations to run
		self.popSize=100 #size of population
		self.mu= 100#individuals selected
		self.lambd= 200#children produced
		self.bitFlipProb=.1 # prob of flipping bits inside mutation
		self.crossoverProb=.2 # prob of crossing over a particular parent
		self.mutationProb=.2 # prob of mutating a particular parent
		self.async=False # run in asynchronous mode
		self.iters=100 #number of simulations to try in asynchronous mode
		self.complexPenalty=False #penalize models with increased complexity
		self.genSteps=10000 # steps to find steady state with fake data
		self.sigmaNetwork=0
		self.sigmaNode=0
		self.hofSize=10
		self.cells=100

def fuzzyAnd(num1, num2, index1, index2,corrMat):
	return min(num1,num2)
def fuzzyOr(num1, num2, index1, index2,corrMat):
	return max(num1, num2)
def naiveAnd(num1, num2, index1, index2,corrMat):
	return num1*num2	
def naiveOr(num1, num2, index1, index2,corrMat):
	return (num1+num2-(num1*num2))
def propOr(A,B,AgivenB, BgivenA):
	return max(0,min(1,A+B-((B*AgivenB+A*BgivenA)/2)))
def propAnd(A,B,AgivenB, BgivenA):
	return max(0,min(1,(A*BgivenA+B*AgivenB)/2))
def Inv(x, inverter): #inverts if necessary then applies hill fun
	if inverter:
		return (1-x)
	else:
		return (x)

class simulatorClass:
	def __init__(self,simTyping):
		self.simType=simTyping
		if simTyping=='prop':
			self.And=propAnd
			self.Or=propOr
			self.corrMat={}
		if simTyping=='fuzzy':
			self.And=fuzzyAnd
			self.Or=fuzzyOr
			self.corrMat=0
		if simTyping=='propNaive':
			self.And=naiveAnd
			self.Or=naiveOr
			self.corrMat=0

def updateNode(currentNode,oldValue,nodeIndividual, model,simulator):
	# we update node by updating shadow and nodes then combining them to update or nodes. 
	andNodes=model.andNodeList[currentNode] # find the list of shadow and nodes we must compute before computing value of current nodes
	andNodeInvertList=model.andNodeInvertList[currentNode] #find list of lists of whether input nodes need to be inverted (corresponds to inputOrder)
	if andLenList[currentNode]==0:
		return oldValue[currentNode] #if no inputs, maintain value
	elif len(andNodes)==1: 
		#if only one input, then can either affect or not affect the node. so either keep the value or update to the single input's value
		if individual[0]==1:
			#if only one input, then set to that number
			value=Inv(oldValue[andNodes[0][0]],andNodeInvertList[0][0])
		else:
			value=oldValue[currentNode] #if no inputs, maintain value
		return value
	else:
		#update nodes with more than one input
		# update and then or
		
		upstreamVals=[]
		for upstream in range(0,len(inputOrder)):
			upstreamVals.append(Inv(oldValue[inputOrder[upstream]],andNodeInvertList[upstream]))
		counter =0
		while counter < len(logicOperatorFlags) and counter+1<len(inputOrder):
			if logicOperatorFlags[counter]==0:
				tempVal=simulator.And(upstreamVals[counter],upstreamVals[counter+1],inputOrder[counter], inputOrder[counter+1],simulator.corrMat)
				inputOrder.pop(counter)
				inputOrder.pop(counter)
				logicOperatorFlags.pop(counter)
				upstreamVals.pop(counter)
				upstreamVals.pop(counter)
				upstreamVals.insert(counter,tempVal)
			else:
				counter=counter+1
			#first one uses the initial logic operator flag to decide and vs or then combines the first two inputs
		while len(upstreamVals)>1:
			tempVal=simulator.Or(upstreamVals.pop(0),upstreamVals.pop(0),2,1,simulator.corrMat)
			upstreamVals.insert(0,tempVal)
				# print(upstreamVals)
		return upstreamVals[0]		

#run a simulation given a starting state
def runModel(individual, params, model, simulator, initValues):
	# do simulation. individual specifies the particular logic rules on the model. params is a generic holder for simulation parameters. 
	
	# set up data storage for simulation, add step 0
	newValue=list(initValues)
	simData=[]
	simData.append(list(newValue))

	# set up the sequence of nodes to be updated
	seq=range(0,len(model.individualParse))
	for node in range(0,len(model.earlyEvalNodes)):
		newValue[model.earlyEvalNodes[node]]=updateNode(model.earlyEvalNodes[node],newValue,individual,individualParse[seq[i]], model,simulator)
	for step in range(0,params.simSteps):
		oldValue=list(newValue)
		if params.async:
			shuffle(seq)
		for i in range(0,len(model.individualParse)):	
			if params.async:
				temp=updateNode(seq[i],newValue,individual, model.individualParse[seq[i]],  model,simulator)
			else:
				temp=updateNode(seq[i],oldValue,individual, model.individualParse[seq[i]], model,simulator)
			newValue[seq[i]]=temp

		simData.append(list(newValue))
	avg= [0 for x in range(0,len(newValue))]
	stdev= [0 for x in range(0,len(newValue))]
	for step in range(0,len(simData)):
		for element in range(0,len(avg)):
			simData[step][element]=simData[step][element]+params.sigmaNetwork*random()
	for step in range(len(simData)-10,len(simData)):
		for element in range(0,len(avg)):
			avg[element]=avg[element]+simData[step][element]
	for element in range(0,len(avg)):
		avg[element]=avg[element]/10
	#for step in range(len(simData)-10,len(simData)):
		#for element in range(0,len(avg)):
			#stdev[element]=stdev[element]+(simData[step][element]-avg[element])^2
	#return (avg, stdev)
	return avg

#run a simulation and average it over iters trials
def averageResultModelSim(individual, params, model, simulator, initValues, iters):
	sum=[0 for x in range(0,len(initValues))]
	for i in range(0,iters):
		avg=runModel( individual, params, model, simulator, initValues,False)
		for j in range(0,len(sum)):
			sum[j]=sum[j]+avg[j]
	avgs=list(sum)
	for i in range(0,len(sum)):
		avgs[i]=sum[i]/float(iters)
	return avgs



