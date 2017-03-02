from gensim.models import Word2Vec
from collections import defaultdict
from time import time
import argparse
import sys
import cPickle

#######Initialize argument parser#######
parser = argparse.ArgumentParser()
parser.add_argument("--model", help="name of word2vec file", default="./data/pubmedpmcwiki")
parser.add_argument("--pairs", help="name of pairs file", default="pairs.txt")
parser.add_argument("--out", help="name of output file", default="output.txt")

def findSimDistance(model, pairsFn, cutoff = .00005):
	simDistDict = defaultdict(list)
	pairsList = readPairs(pairsFn, model)
	for index, newWord in enumerate(model.index2word):
		#sys.stdout.write("Percent done: %.2f%%\r" % float(index/ len(model.index2word) * 100))
		sys.stdout.write("Percent done: " + str((float(index)/len(model.index2word))*100.0)+"%%\r")
		sys.stdout.flush()
		for pair in pairsList:
			word1 = pair[0]
			word2 = pair[1]
			dist = pair[2]
			newDist = model.similarity(word1, newWord)
			if abs(newDist - dist) <= cutoff and word1 != newWord and word2 != newWord : simDistDict[pair].append((newWord, newDist))
	with open("simDistDict", 'wb') as f:
		cPickle.dump(simDistDict, f, cPickle.HIGHEST_PROTOCOL)
	return simDistDict

def readPairs(pairsFn, model):
	pairsList = []
	with open(pairsFn, 'r') as pairs:
		for line in pairs:
			word1 = line.split(" ")[0].rstrip('\n')
			word2 = line.split(" ")[1].rstrip('\n')
			dist = model.similarity(word1, word2)
			pairsList.append((word1, word2, dist))
	return pairsList

def computeAll(distDict, modelFn, pairsFn, outputFn):
	with open(outputFn, 'a') as out:
		out.write("This file was generated using model " + modelFn + " and pair file " + pairsFn + ".\n")
		out.write("The Cosine Distance cutoff for this file was 0.00005.\n")
		for pairTuple, newWordTupleList in distDict.items():
			out.write(str(pairTuple[0]) + " and " + str(pairTuple[1]) + ": " + str(round(pairTuple[2],5)) + "\n")
			out.write("----------------\n")
			for newWord, newDist in newWordTupleList:
				out.write("\t" + str(newWord) + ": " + str(round(float(newDist),5)) + "\n")

if __name__ == '__main__':
	flags = parser.parse_args()
	modelFn = flags.model
	pairsFn = flags.pairs
	outFn = flags.out
	print("Using word2vec model " + modelFn)
	print("Reading from " + pairsFn)
	print("Writing to " + outFn)
	print("Loading model.....")
	start = time()
	model = Word2Vec.load(modelFn, mmap = 'r')
	end = time()-start
	print("Done loading model! It took " + str(round(end,2)) + " seconds!")
	print("Computing distance maps for all pairs...")
	start = time()
	simDistDict = findSimDistance(model, pairsFn)
	end = time() - start
	print("Finished computing distance maps! It took " + str(round(end, 2)) + " seconds!")
	#simDistDict = cPickle.load(open("simDistDict"))
	print("Writing to output file....")
	computeAll(simDistDict, modelFn, pairsFn, outFn)
	print("Done writing to output file!")
	
	
