from gensim.models import Word2Vec
from collections import defaultdict
from time import time
import argparse
import sys
import cPickle
from pandas import *
import numpy as np
import random
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import normalize
from sklearn import decomposition
from sklearn.neural_network import MLPClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import VotingClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib
#import matplotlib.pyplot as plt

#######Initialize argument parser#######
parser = argparse.ArgumentParser()
parser.add_argument("--model", help="name of word2vec file", default="./data/pubmedpmcwiki")
parser.add_argument("--infn", help="name of interaction file", default="interactions.tsv")
parser.add_argument("--outfn", help="name of output file", default="outputNicoleTop.txt")
parser.add_argument("--diff", help="if True take the difference of drug and gene word vectors", default=False)
parser.add_argument("--aim1", default = "Nicole_Predictions.tsv")
parser.add_argument("--plot", default = True)
def readInteractionFile(interactionsFn):
	interactionDict = defaultdict(list)
	geneSet = set([])
	drugSet = set([])
	interactions = read_table(interactionsFn)
	for row in range(len(interactions)):
		geneName = interactions['entrez_gene_symbol'][row]
		drugName = interactions['drug_primary_name'][row]
		interactionDict[geneName].append(drugName)
		geneSet.add(geneName)
		drugSet.add(drugName)
	return (interactionDict, geneSet, drugSet)

def buildNegativeSamples(interactionDict, geneSet, drugSet, numNegSamples, model, diff):
	random.seed(2443)
	numGenes = len(geneSet)
	negSamples = list()
	for i in range(numNegSamples):
		while True:
			geneName = random.sample(geneSet, 1)[0]
			drugName = random.sample(drugSet, 1)[0]
			while (drugName in interactionDict[geneName]): drugName = random.sample(drugSet, 1)
			try:
				geneVector = model[geneName]
				drugVector = model[drugName]
				break
			except:
				continue
		geneVec = np.array(geneVector)
		drugVec = np.array(drugVector)
		difference = np.array(drugVec - geneVec)
		if len(difference) != 200: 
			i -= 1
			continue
		if diff: negSamples.append(difference)
		else: negSamples.append(np.append(geneVec, drugVec))
	return negSamples

def buildPosSamples(interactionDict, model, diff):
	posSamples = list()
	posDrugs = list()
	posGenes = list()
	posVec = list()
	for geneName, drugList in interactionDict.items():
		for drugName in drugList:
			try:
				geneVector = model[geneName]
				drugVector = model[drugName]
			except:
				continue
			geneVec = np.array(geneVector)
			drugVec = np.array(drugVector)
			difference = np.array(drugVec - geneVec)
			if diff: posSamples.append(difference)
			else: posSamples.append(np.append(geneVec, drugVec))

	return posSamples

def mergePosNeg(posSamples, negSamples):
	with open("posSamples", 'wb') as f:
		cPickle.dump(posSamples, f, cPickle.HIGHEST_PROTOCOL)
	with open("negSamples", 'wb') as f:
		cPickle.dump(negSamples, f, cPickle.HIGHEST_PROTOCOL)
	totalSamples = np.vstack((posSamples, negSamples))
	totalLabels = np.vstack((np.ones((len(posSamples), 1)),np.zeros((len(negSamples), 1))))
	with open("totalSamples", 'wb') as f:
		cPickle.dump(totalSamples, f, cPickle.HIGHEST_PROTOCOL)
	with open("totalLabels", 'wb') as f:
		cPickle.dump(totalLabels, f, cPickle.HIGHEST_PROTOCOL)
	totalSamples = normalize(totalSamples, axis=0)
	c, r = totalLabels.shape
	totalLabels = totalLabels.reshape(c,)
	return totalSamples, totalLabels

def crossValSVM(totalSamples, totalLabels):
	# pca = decomposition.PCA(n_components = 100)
	# pca.fit(totalSamples)
	# totalSamples = pca.transform(totalSamples)
	# print("Applied PCA to samples. We now have " + str(len(totalSamples[0])) + " features.")
	X_train, X_test, y_train, y_test = train_test_split(totalSamples, totalLabels, test_size=0.25, random_state=2243)
	maxMean = 0
	maxclf = None
	for kernel in ('linear', 'poly', 'rbf'):
		clf = svm.SVC(kernel=kernel, gamma=2)
		scores = cross_val_score(clf, totalSamples, totalLabels, cv=5)
		print("Accuracy for " + kernel + " kernel SVM: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
		if scores.mean() > maxMean:
			maxMean = scores.mean()
			maxclf = clf
	return clf, maxMean

def crossValMLP(totalSamples, totalLabels):
	clf = MLPClassifier(solver='lbfgs', alpha = 1e-5, hidden_layer_sizes=(5,2), random_state = 2443)
	scores = cross_val_score(clf, totalSamples, totalLabels, cv=5)
	print("Accuracy for mlp-classifier: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
	return clf, scores.mean()	

def crossValKNN(totalSamples, totalLabels, n):
	clf = KNeighborsClassifier(n_neighbors=n)
	scores = cross_val_score(clf, totalSamples, totalLabels, cv=5)
	print("Accuracy for knn: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
	return clf, scores.mean()

def crossValDecisionTree(totalSamples, totalLabels):
	clf = DecisionTreeClassifier(max_depth=4)
	scores = cross_val_score(clf, totalSamples, totalLabels, cv=5)
	print("Accuracy for decision tree: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
	return clf, scores.mean()

def crossValRandomForest(totalSamples, totalLabels):
	clf = RandomForestClassifier(n_estimators=10, max_depth=None, min_samples_split=2, random_state=0)
	scores = cross_val_score(clf, totalSamples, totalLabels, cv=5)
	print("Accuracy for random forest: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
	return clf, scores.mean()

def crossValVoting(totalSamples, totalLabels, estimators, votingWeights):
	clf = VotingClassifier(estimators=estimators, voting = 'hard')
	scores = cross_val_score(clf, totalSamples, totalLabels, cv=5)
	print("Accuracy for voting: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
	return clf, scores.mean()

def crossValBoost(totalSamples, totalLabels):
	clf = AdaBoostClassifier(n_estimators=100)
	scores = cross_val_score(clf, totalSamples, totalLabels, cv=5)
	print("Accuracy for boosting: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
	return clf, scores.mean()

def evalPredictors():
	flags = parser.parse_args()
	modelFn = flags.model
	interactionsFn = flags.infn
	outFn = flags.outfn
	diff = flags.diff
	print("Using word2vec model " + modelFn)
	print("Reading from " + interactionsFn)
	print("Writing to " + outFn)
	if diff: print ("Taking the difference of the drug and gene word vectors")
	print("Loading model.....")
	start = time()
	model = Word2Vec.load(modelFn, mmap = 'r')
	end = time()-start
	print("Done loading model! It took " + str(round(end,2)) + " seconds!")
	interactionDict, geneSet, drugSet = readInteractionFile(interactionsFn)
	print("Building positive samples...")
	pos = buildPosSamples(interactionDict, model, diff)
	print("Done building positive samples!")
	print("Building negative samples...")
	neg = buildNegativeSamples(interactionDict, geneSet, drugSet, len(pos), model, diff)
	print("Done building negative samples!")
	totalSamples, totalLabels = mergePosNeg(pos, neg)
	print("We now have %d total samples"%len(totalSamples))
	svm, svmMean = crossValSVM(totalSamples, totalLabels)
	mlp, mlpMean = crossValMLP(totalSamples, totalLabels)
	knn, knnMean = crossValKNN(totalSamples, totalLabels, 25)
	dt, dtMean = crossValDecisionTree(totalSamples, totalLabels)
	rf, rfMean = crossValRandomForest(totalSamples, totalLabels)
	ab, abMean = crossValBoost(totalSamples, totalLabels)
	modelList = [mlpMean, knnMean, dtMean, rfMean, abMean]
	weights = [int(x/max(modelList)) + 1 for x in modelList]
	estimators = [("svm", svm), ("mlp",mlp), ("knn",knn), ("dt", dt), ("rf", rf), ("ab", ab)]
	weights = [x/x for x in modelList]
	voting, votingMean = crossValVoting(totalSamples, totalLabels, estimators, weights)

def evalAim1():
	flags = parser.parse_args()
	aim1Fn = flags.aim1
	outFn = flags.outfn
	modelFn = flags.model
	interactionsFn = flags.infn
	diff = flags.diff
	print("Evaluating predictions in %s"%aim1Fn)
	print("Writing output to %s"%outFn)
	print("Loading model.....")
	start = time()
	model = Word2Vec.load(modelFn, mmap = 'r')
	end = time()-start
	print("Done loading model! It took " + str(round(end,2)) + " seconds!")
	interactionDict, geneSet, drugSet = readInteractionFile(interactionsFn)
	print("Building positive samples...")
	pos = buildPosSamples(interactionDict, model, diff)
	print("Done building positive samples!")
	print("Building negative samples...")
	neg = buildNegativeSamples(interactionDict, geneSet, drugSet, len(pos), model, diff)
	print("Done building negative samples!")
	totalSamples, totalLabels = mergePosNeg(pos, neg)
	print("We have %d total samples!"%len(totalSamples))


	data = pandas.read_table(aim1Fn)
	drugList = data['Drug']
	geneLists = data['Gene']

	clf = KNeighborsClassifier(n_neighbors=25)
	clf.fit(totalSamples, totalLabels)
	print("We have %d total samples!"%len(totalSamples))
	predictions = defaultdict(float)

	print("Generating predictions...")
	with open(outFn, 'w') as out:
		for index, drug in enumerate(drugList):
			geneList = geneLists[index].split(",")
			numGenes = 0
			for gene in geneList:
				try:
					features = np.append(np.array(model[gene]), np.array(model[drug]))

				except:
					print gene
					print drug
					continue
				prediction = clf.predict(features.reshape(1, -1))
				predictions[drug] += prediction
				numGenes += 1
			
			if numGenes > 0: out.write("%s and %s: %f\n"%(drug, gene, predictions[drug]/numGenes))


	print("Done generating predictions!")

if __name__ == '__main__':
	evalAim1()



