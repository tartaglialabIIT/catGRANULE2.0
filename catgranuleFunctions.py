import numpy as np
from matplotlib import pyplot as plt
from Bio.Seq import Seq
import ast
from numpy import nan

from glob import glob

import math
from sklearn.decomposition import PCA

# from dtw import dtw
from scipy.signal import savgol_filter

import subprocess
import os
import pandas as pd
import random
from Bio.Seq import Seq
from Bio import SeqIO
from itertools import cycle

from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import roc_auc_score

# Import some data to play with
from scipy.stats import pearsonr
import requests as r
from io import StringIO

from collections import Counter
import argparse
import sys
#from stringScalesFunctions import *
import csv
import joblib
import re
import time
import gzip

from Bio import PDB
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB import Selection

#import pymol2
import multiprocessing as mp

from matplotlib import pyplot as plt
import seaborn as sns
import math
from sklearn.decomposition import PCA
import pickle
import random
#from dtw import dtw
from scipy.signal import savgol_filter

import random

from itertools import cycle

from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_auc_score

# Import some data to play with
from scipy.stats import pearsonr

from collections import Counter
from mpl_toolkits.mplot3d import Axes3D

from sklearn.metrics import accuracy_score
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
import requests

from matplotlib.colors import Normalize

from Bio.PDB.DSSP import DSSP
from Bio.PDB import is_aa

import json



# HumanProt = list(SeqIO.parse("HumanProteome.fasta", "fasta"))

def chunkstring(string, length):
    return [string[0 + i:length + i] for i in range(0, len(string), length)]


def Stringentropy(string):
    "Calculates the Shannon entropy of a string"

    # get probability of chars in string
    prob = [float(string.count(c)) / len(string) for c in dict.fromkeys(list(string))]

    # calculate the entropy
    entropy = - sum([p * math.log(p) / math.log(2.0) for p in prob])

    return entropy


def convert(s):
    # initialization of string to ""
    new = ""

    # traverse in the string
    for x in s:
        new += x

        # return string
    return new


def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def smooth2(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='valid')
    return y_smooth


def remove_values_from_list(the_list, val):
    return [value for value in the_list if value != val]


def condensateMixtList(expT):
    return [np.concatenate(expT[i], axis=0) for i in range(np.size(expT))]


manhattan_distance = lambda x, y: np.abs(x - y)


def computeDeltaChemPhysPCA(ScaleNames, PCAcomp, WTseq, mutSeq):
    first = PCAcomp
    count = 0
    newScaleWT = []

    newScaleMut = []

    for fn in ScaleNames:
        newScaleWT.append(first[count] * np.array(transformSeq(fn, WTseq)))

        newScaleMut.append(first[count] * np.array(transformSeq(fn, mutSeq)))

        count = count + 1
    amin_1 = np.array(newScaleWT)

    amin_2 = np.array(newScaleMut)

    return np.mean(amin_1 - amin_2)


def computeMutationScorePCA(ScaleNames, PCAcomp, WTseq, mutation):
    numberMutations = np.shape(mutation)[0]
    mutScore = []
    bindingScore = []
    for i in range(numberMutations):
        # print(i)

        index, wt, mut = mutation.iat[i, 0], mutation.iat[i, 1], mutation.iat[i, 2]
        bindingScore.append(mutation.iat[i, 3])
        mutSeq = seqMutation(WTseq, index, str(mut))

        deltaMut = computeDeltaChemPhysPCA(ScaleNames, PCAcomp, WTseq, mutSeq)
        mutScore.append(deltaMut)

    return bindingScore, mutScore


def transformSeq(filename, seq):
    fileT = open('seqTemporanea.txt', 'w')
    fileT.write(seq)
    #     np.savetxt('seqTemporanea.txt',seq)
    fileT.close()

    tocall = 'bash ' + filename + ' seqTemporanea.txt -> TemporaneaChemPhys.txt'
    os.system(tocall)

    # fileT2= open('TemporaneaChemPhys.txt','r')

    # ChemPhys_Seq = [i for i in fileT2.readlines()]

    ChemPhys_Seq = np.loadtxt('TemporaneaChemPhys.txt')
    os.system('rm seqTemporanea.txt')

    os.system('rm TemporaneaChemPhys.txt')

    return ChemPhys_Seq


def computeDeltaChemPhys(filename, WTseq, mutSeq):
    amin_1 = transformSeq(filename, WTseq)[:, 1]

    amin_2 = transformSeq(filename, mutSeq)[:, 1]

    return np.mean(amin_1 - amin_2)


def seqMutation(seq, mutationIndex, mutatedAmino):
    mutSeq = seq[:mutationIndex] + mutatedAmino + seq[mutationIndex + 1:]
    return mutSeq


def computeMutationScore(filename, WTseq, mutation):
    fn = filename
    numberMutations = np.shape(mutation)[0]
    mutScore = []
    bindingScore = []
    for i in range(numberMutations):
        # print(i)

        index, wt, mut = mutation.iat[i, 0], mutation.iat[i, 1], mutation.iat[i, 2]
        bindingScore.append(mutation.iat[i, 3])
        mutSeq = seqMutation(WTseq, index, str(mut))

        deltaMut = computeDeltaChemPhys(fn, WTseq, mutSeq)
        mutScore.append(deltaMut)

    return bindingScore, mutScore


def giveIndx(mutScore, tresh):
    return [(np.sign(i - tresh) + 1) / 2 for i in mutScore]


def diff_letters(a, b):
    return sum(a[i] != b[i] for i in range(len(a)))


def translate_from_dict(original_text, dictionary_of_translations):
    out = original_text
    for target in dictionary_of_translations:
        trans = str.maketrans(target, dictionary_of_translations[target])
        out = out.translate(trans)
    return out


def flipIndx(arr):
    return [1 if i == 0 else 0 for i in arr]


aminoList = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def RandomMix(seq, perc):
    seq = chunkstring(seq, 1)
    seqNew = []
    for i in seq:
        rn = np.random.random()
        if perc > rn:
            j = aminoList[np.random.randint(len(aminoList))]
        else:
            j = i

        seqNew.append(j)

    ns = ''.join(seqNew)
    return ns


def multipleRandomSeq(seq, N, perc):
    seqList = []
    for i in range(N):
        seqList.append(RandomMix(seq, perc))

    return seqList


def multipleRandomSeq_transform(fn, seq, N, perc):
    seqList = multipleRandomSeq(seq, N, perc)
    scaleList = []
    for j in range(len(seqList)):
        scaledS = transformSeq(fn, seqList[j])
        scaleList.append(scaledS)
    return scaleList


def multipleRandomSeq_transformMeanScore(fn, seq, N, perc):
    seqList = multipleRandomSeq(seq, N, perc)
    scaleList = []
    for j in range(len(seqList)):
        scaledS = transformSeq(fn, seqList[j])
        scaleList.append(np.mean(scaledS))
    return scaleList


def TopBotRND(aggScore, scoreList, fn, WtSeq_utile, perc, plot=False):
    sortAgg = np.sort(aggScore)
    indx_sort = np.argsort(aggScore)
    scorSort = np.array(scoreList)[indx_sort]

    topBotList = [10, 200, 500, 2000, 5000, int(np.size(sortAgg) / 2)]

    topBotList = [10, 30, 50, 70, 100, 150, 200, 300, 400, 500, int(np.size(sortAgg) / 2)]

    roc_auc = []

    # perc =0.05
    if plot:
        f, ax = plt.subplots(1, 2, figsize=(30, 30))

    for i in topBotList:

        top = sortAgg[:i]
        bot = multipleRandomSeq_transformMeanScore(fn, WtSeq_utile, i, perc)

        topScore = scorSort[:i]
        botScore = np.zeros(len(scorSort[-i:])) - 1

        newAgg = np.concatenate([np.array(top), np.array(bot)])
        newScore = np.concatenate([np.array(topScore), np.array(botScore)])

        indx_agg = flipIndx(giveIndx(newScore, 0))

        x, y, q = roc_curve(indx_agg, newAgg)
        auc_val = roc_auc_score(indx_agg, newAgg)
        roc_auc.append(auc_val)

        if plot:
            ax[0].plot(x, y, label='AUC = %0.2f  top-bottom = %.4f' % (auc_val, i))
            ax[0].legend(loc='best')

    if plot:
        ax[1].plot(topBotList, roc_auc)
        ax[1].set_xlabel('Top-Bottom size')
        ax[1].set_ylabel('AUC')
        plt.show()

    return roc_auc, topBotList


# indx on experimental
def TopBot(aggScore, scoreList, fn, WtSeq_utile, plot=False):
    '''
    indx onexperimental
    top and bottom on theory
    '''

    sortAgg = np.sort(aggScore)
    indx_sort = np.argsort(aggScore)
    scorSort = np.array(scoreList)[indx_sort]

    topBotList = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 500, 600, 750, 1000, 1500, 2000,
                  3000, 4500, 5000, 7500, 10000, 12000, 15000, int(np.size(sortAgg) / 2)]

    topBotList = [10, 30, 50, 70, 100, 150, 200, 300, 400, 500, int(np.size(sortAgg) / 2)]

    roc_auc = []
    if plot:
        f, ax = plt.subplots(1, 2, figsize=(30, 30))

    for i in topBotList:

        top = sortAgg[:i]
        bot = sortAgg[-i:]

        topScore = scorSort[:i]
        botScore = scorSort[-i:]

        newAgg = np.concatenate([np.array(top), np.array(bot)])
        newScore = np.concatenate([np.array(topScore), np.array(botScore)])

        indx_agg = flipIndx(giveIndx(newScore, 0))

        x, y, q = roc_curve(indx_agg, newAgg)
        auc_val = roc_auc_score(indx_agg, newAgg)
        roc_auc.append(auc_val)

        if plot:
            ax[0].plot(x, y, label='AUC = %0.2f  top-bottom = %.4f' % (auc_val, i))
            ax[0].legend(loc='best')

    if plot:
        ax[1].plot(topBotList, roc_auc, '.')

        ax[1].plot(topBotList, roc_auc, '-')

        ax[1].set_xlabel('Top-Bottom size')
        ax[1].set_ylabel('AUC')
        ax[1].set_xscale('log')

    return roc_auc, topBotList


# indx on experimental
def TopBot2(aggScore, scoreList, fn, WtSeq_utile, plot=False):
    '''
    indx on experimental
    top and bottom on experimental
    '''

    sortAgg = np.sort(aggScore)
    indx_sort = np.argsort(aggScore)
    scorSort = np.array(scoreList)[indx_sort]

    scorSort = np.sort(scoreList)
    indx_sort = np.argsort(scoreList)
    sortAgg = np.array(sortAgg)[indx_sort]

    topBotList = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 500, 600, 750, 1000, 1500, 2000,
                  3000, 4500, 5000, 7500, 10000, 12000, 15000, int(np.size(sortAgg) / 2)]

    # topBotList =[100,500,1000,5000,15000, int(np.size(sortAgg)/2)]

    roc_auc = []
    if plot:
        f, ax = plt.subplots(1, 2, figsize=(30, 30))

    for i in topBotList:

        top = sortAgg[:i]
        bot = sortAgg[-i:]

        topScore = scorSort[:i]
        botScore = scorSort[-i:]

        newAgg = np.concatenate([np.array(top), np.array(bot)])
        newScore = np.concatenate([np.array(topScore), np.array(botScore)])

        indx_agg = flipIndx(giveIndx(newScore, 0))

        x, y, q = roc_curve(indx_agg, newAgg)
        auc_val = roc_auc_score(indx_agg, newAgg)
        roc_auc.append(auc_val)

        if plot:
            ax[0].plot(x, y, label='AUC = %0.2f  top-bottom = %.4f' % (auc_val, i))
            ax[0].legend(loc='best')

    if plot:
        ax[1].plot(topBotList, roc_auc, '.')

        ax[1].plot(topBotList, roc_auc, '-')

        ax[1].set_xlabel('Top-Bottom size')
        ax[1].set_ylabel('AUC')
        ax[1].set_xscale('log')

    return roc_auc, topBotList


def dataVsRND(aggScore, sizeList, fn, WtSeq_utile, perc, plot=False):
    roc_auc = []

    # perc =0.05
    if plot:
        f, ax = plt.subplots(1, 2, figsize=(30, 30))

    for i in sizeList:

        sortAgg = random.sample(aggScore, i)
        top = sortAgg
        bot = multipleRandomSeq_transformMeanScore(fn, WtSeq_utile, i, perc)

        topScore = np.zeros(len(sortAgg)) + 1
        botScore = np.zeros(len(sortAgg))

        newAgg = np.concatenate([np.array(top), np.array(bot)])
        indx_agg = np.concatenate([np.array(topScore), np.array(botScore)])

        # indx_agg = flipIndx(giveIndx(newScore,0))

        x, y, q = roc_curve(indx_agg, newAgg)
        auc_val = roc_auc_score(indx_agg, newAgg)
        roc_auc.append(auc_val)

        if plot:
            ax[0].plot(x, y, label='AUC = %0.2f  Pool size = %.4f' % (auc_val, i))
            ax[0].legend(loc='best')

    if plot:
        ax[1].plot(sizeList, roc_auc)
        ax[1].set_xlabel('Pool size')
        ax[1].set_ylabel('AUC')
        plt.show()

    return roc_auc, sizeList


def RandomMix_num(seq, num):
    seq = chunkstring(seq, 1)
    seqNew = []
    for i in range(num):
        rn = np.random.randint(len(seq))

        seq[rn] = aminoList[np.random.randint(len(aminoList))]

    ns = ''.join(seq)
    return ns


def RandomMix_num_2(seq, num):
    seq = chunkstring(seq, 1)
    seqNew = []
    newAmList = []
    oldAmList = []
    positionMutation = []
    for i in range(num):
        rn = np.random.randint(len(seq))

        newAm = aminoList[np.random.randint(len(aminoList))]
        oldAmList.append(seq[rn])
        newAmList.append(newAm)
        positionMutation.append(rn)

        seq[rn] = newAm

    ns = ''.join(seq)
    return ns, oldAmList, newAmList, positionMutation


def multipleRandomSeq_num(seq, N, num):
    seqList = []
    for i in range(N):
        seqList.append(RandomMix_num(seq, num))

    return seqList


def multipleRandomSeq_transform_num(fn, seq, N, num):
    seqList = multipleRandomSeq_num(seq, N, num)
    scaleList = []
    for j in range(len(seqList)):
        scaledS = transformSeq(fn, seqList[j])
        scaleList.append(scaledS)
    return scaleList


def multipleRandomSeq_transformMeanScore_num(fn, seq, N, num):
    seqList = multipleRandomSeq_num(seq, N, num)
    scaleList = []
    for j in range(len(seqList)):
        scaledS = transformSeq(fn, seqList[j])
        scaleList.append(np.mean(scaledS))
    return scaleList


def dataVsRND_num(aggScore, mutList, indxList, fn, WtSeq_utile, plot=False):
    roc_auc = []

    # perc =0.05
    if plot:
        f, ax = plt.subplots(1, 2, figsize=(30, 30))

    for i in range(len(mutList)):
        print(len(indxList[i]))
        if (len(indxList[i])) > 0:
            sortAgg = aggScore[indxList[i]]
            top = sortAgg
            bot = multipleRandomSeq_transformMeanScore_num(fn, WtSeq_utile, len(sortAgg), i)

            topScore = np.zeros(len(sortAgg)) + 1
            botScore = np.zeros(len(sortAgg))

            newAgg = np.concatenate([np.array(top), np.array(bot)])
            indx_agg = np.concatenate([np.array(topScore), np.array(botScore)])

            # indx_agg = flipIndx(giveIndx(newScore,0))

            x, y, q = roc_curve(indx_agg, newAgg)
            auc_val = roc_auc_score(indx_agg, newAgg)
            roc_auc.append(auc_val)

        if plot:
            ax[0].plot(x, y, label='AUC = %0.2f  N mutations = %.4f' % (auc_val, i))
            ax[0].legend(loc='best')

    if plot:
        ax[1].plot(range(len(mutList)), roc_auc)
        ax[1].set_xlabel('Pool size')
        ax[1].set_ylabel('AUC')
        plt.show()

    return roc_auc, range(len(mutList))


def dataVsRND_num_poolsize(aggScore, mutList, indxList, poolsize, fn, WtSeq_utile, plot=False):
    roc_auc = []

    # perc =0.05
    if plot:
        f, ax = plt.subplots(1, 2, figsize=(30, 30))
    True_IL = []

    for i in range(len(mutList)):
        if i % 100 == 0:
            print(i / len(mutList))
        if (len(indxList[i])) > 0:
            # print('NumberMutations:  ',i)

            sortAgg = aggScore[indxList[i]]
            if len(sortAgg) > poolsize:
                sortAgg = random.sample(list(sortAgg), poolsize)

            top = sortAgg
            bot = multipleRandomSeq_transformMeanScore_num(fn, WtSeq_utile, len(sortAgg), i)

            topScore = np.zeros(len(sortAgg)) + 1
            botScore = np.zeros(len(sortAgg))

            newAgg = np.concatenate([np.array(top), np.array(bot)])
            indx_agg = np.concatenate([np.array(topScore), np.array(botScore)])

            # indx_agg = flipIndx(giveIndx(newScore,0))

            x, y, q = roc_curve(indx_agg, newAgg)
            auc_val = roc_auc_score(indx_agg, newAgg)
            roc_auc.append(auc_val)
            True_IL.append(i)
        if plot:
            ax[0].plot(x, y, label='AUC = %0.2f  N mutations = %.d' % (auc_val, i))
            # ax[0].legend(loc='best')

    if plot:
        ax[1].plot(True_IL, roc_auc)
        ax[1].set_xlabel('Number of mutations')
        ax[1].set_ylabel('AUC')
        f.suptitle('Pool-Size %d' % poolsize)
        plt.show()

    return roc_auc, True_IL


def dataVsRND_num_poolsize0(aggScore, mutList, indxList, poolsize, fn, WtSeq_utile, plot=False):
    roc_auc = []

    # perc =0.05
    if plot:
        f, ax = plt.subplots(1, 2, figsize=(30, 30))
    True_IL = []

    for i in range(len(mutList)):
        if i % 100 == 0:
            print(i / len(mutList))

        # print('NumberMutations:  ',i)

        sortAgg = aggScore[indxList[i]]
        if len(sortAgg) > poolsize:
            sortAgg = random.sample(list(sortAgg), poolsize)

        top = sortAgg
        bot = multipleRandomSeq_transformMeanScore_num(fn, WtSeq_utile, len(sortAgg), i)

        topScore = np.zeros(len(sortAgg)) + 1
        botScore = np.zeros(len(sortAgg))

        newAgg = np.concatenate([np.array(top), np.array(bot)])
        indx_agg = np.concatenate([np.array(topScore), np.array(botScore)])

        # indx_agg = flipIndx(giveIndx(newScore,0))

        x, y, q = roc_curve(indx_agg, newAgg)
        auc_val = roc_auc_score(indx_agg, newAgg)
        roc_auc.append(auc_val)
        True_IL.append(i)
        if plot:
            ax[0].plot(x, y, label='AUC = %0.2f  N mutations = %.d' % (auc_val, i))
            # ax[0].legend(loc='best')

    if plot:
        ax[1].plot(True_IL, roc_auc)
        ax[1].set_xlabel('Number of mutations')
        ax[1].set_ylabel('AUC')
        f.suptitle('Pool-Size %d' % poolsize)
        plt.show()

    return roc_auc, True_IL


def mutationDistribution(filename, WtSeq_utile):
    f = open(filename, 'r')
    lines = f.readlines()

    difflet = []

    for i in lines:
        mutSq = list(i[:-1])

        if len(mutSq) == len(WtSeq_utile):
            nl = diff_letters(WtSeq_utile, mutSq)
            # if nl<15:
            #   difflet.append(nl)

            difflet.append(nl)

    un = np.unique(difflet)

    mutList = []
    indxList = []
    for i in range(np.max(un)):
        mutList.append([])
        indxList.append([])

    count = 0
    for i in lines:
        mutSq = list(i[:-1])

        if len(mutSq) == len(WtSeq_utile):
            nl = diff_letters(WtSeq_utile, mutSq)

            mutList[nl - 1].append(mutSq)
            indxList[nl - 1].append(count)

            count = count + 1

    mutSize = []
    for i in mutList:
        mutSize.append(len(i))

    return mutSize / np.sum(mutSize), mutList, indxList


# compute Right AggScore from sequences with same lenght

def aggSameLenght(aggScore, WtSeq_utile):
    aggScore2 = []
    aggScore3 = []
    for i in range(len(aggScore)):
        if len(WtSeq_utile) - 1 < seqSize[i] <= len(WtSeq_utile):
            aggScore2.append(aggScore[i])

        else:
            aggScore3.append(aggScore[i])

    print(len(aggScore), len(aggScore2), len(aggScore3))

    return aggScore2, aggScore3


def cutTresh(distr, tresh):
    return [1 if i > tresh else 0 for i in distr]


# copmute agg score from fn scale and seqList list of sequences
def computeAggScore(seqList, fn):
    aggList = []

    aggScore = []

    for j in range(np.size(seqList)):

        if j % 1000 == 0:
            print(j / np.size(seqList))

        scaledS = transformSeq(fn, seqList[j])
        aggList.append(scaledS)
        aggScore.append(np.mean(scaledS))

    return aggScore, aggList


def sequenceList(filename):
    f = open(filename, 'r')
    lines = f.readlines()

    # generate the sequence list

    seqList = []
    SeqId = []
    count = 0
    for i in range(len(lines)):
        if i % 2 == 1:
            seqList.append("".join(lines[i][:-1]))
        if i % 2 == 0:
            SeqId.append(lines[i])
        count = count + 1

    f.close()
    print(count)

    return seqList, SeqId


def LenghtDistribution(seqList, plot=False, density=True):
    seqSize = []
    for i in seqList:
        seqSize.append(len(i))

    hist, bins2 = np.histogram(seqSize, bins=50, density=True)

    if plot:
        plt.plot(bins2[:-1], hist, '-')

        plt.xlabel('Sequence Lenght')
        plt.ylabel('Count')
        plt.show()

    return hist, bins2


def ScoreSeq(fn, seq):
    Amin0 = transformSeq(fn, seq)[:, 1]

    return np.mean(Amin0)


def div_d(my_dict):
    sum_p = sum(my_dict.values())

    for i in my_dict:
        my_dict[i] = float(my_dict[i] / sum_p)

    return my_dict


# aminoFreq = []
# oldDict=Counter({})
# for prot in HumanProt[:]:
#     #aminoFreq.append([])
#     seq=str(prot.seq)
#     aminoFreq.append({i:seq.count(i) for i in seq})

#     newDict = Counter({i:seq.count(i) for i in seq})

#     oldDict = oldDict+ newDict

# aminoFreq=div_d(oldDict)
# namesAminoOrdered= list(np.array(list(aminoFreq.items()))[:,0])
# aminoDist=list(aminoFreq.values())


def linearRandom(aminoDist, namesAminoOrdered):
    z = np.random.random()
    pre = 0
    post = 0
    for i in range(0, len(aminoDist)):
        post = post + aminoDist[i]

        if pre < z <= post:
            return namesAminoOrdered[i]
        pre = pre + aminoDist[i]


def RandomMix_num_weighted(seq, num):
    seq = chunkstring(seq, 1)
    seqNew = []
    for i in range(num):
        rn = np.random.randint(len(seq))

        seq[rn] = linearRandom(aminoDist, namesAminoOrdered)

    ns = ''.join(seq)
    return ns


def camFOLD(seq):
    path0 = os.getcwd()
    path = '/Users/mmonti/Documents/Scienza/Mutations_and_functions/gianAlgorithms/camZFOLD/'

    os.chdir(path)

    tocall1 = 'sh start.sh ' + seq + ' > results.txt'

    # tocall1 = 'cd ' + path + ' ; ' + tocall1
    os.system(tocall1)
    f = open('results.txt', 'r')

    r = f.readlines()

    if r[2][44] == '[' or r[2][43] == '[':

        r1 = float(r[2][38:43])

        r2 = float(r[3][38:43])
    else:

        r1 = float(r[2][38:44])

        r2 = float(r[3][38:44])

    f.close()
    os.chdir(path0)

    return r1, r2


def camFOLD_profile(seq):
    path0 = os.getcwd()
    path = '/Users/mmonti/Documents/Scienza/Mutations_and_functions/gianAlgorithms/camZFOLD/'

    os.chdir(path)

    tocall1 = 'sh start.sh ' + seq + ' > results.txt'

    # tocall1 = 'cd ' + path + ' ; ' + tocall1
    os.system(tocall1)
    f = open('results.txt', 'r')

    r = f.readlines()

    x_indx = []

    x_score = []

    for i in r[6:]:
        y = tuple(float(s) for s in i.strip().split(" "))
        x_indx.append(y[0])

        x_score.append(y[1])

    f.close()

    os.chdir(path0)

    return x_indx, x_score


def Ziggregator(seq):
    path0 = os.getcwd()

    path = '/Users/mmonti/Documents/Scienza/Mutations_and_functions/gianAlgorithms/camZAGG/'

    os.chdir(path)

    tocall1 = 'echo ' + seq + ' > sequence.txt; echo "7" >> sequence.txt'
    # tocall1 = 'cd ' + path + ' ; ' + tocall1
    os.system(tocall1)

    tocall2 = ' ./zagg'
    os.system(tocall2)
    tocall3 = "awk '($2>0){s=s+$2;k++} END{print s/k}' zagg7.txt > res.txt"
    os.system(tocall3)

    res = np.loadtxt('res.txt')

    os.system('rm res.txt')

    tocall4 = "awk '{s=s+$2;k++} END{print s/k}' zagg7.txt  > res2.txt"
    os.system(tocall4)

    res2 = np.loadtxt('res2.txt')

    os.system('rm res2.txt')

    os.chdir(path0)

    if res.size == 0:
        res = 0
    if res2.size == 0:
        res2 = 0

    return float(res), float(res2)


def Ziggregator_profile(seq):
    path0 = os.getcwd()

    path = '/Users/mmonti/Documents/Scienza/Mutations_and_functions/gianAlgorithms/camZAGG/'

    os.chdir(path)

    tocall1 = 'echo ' + seq + ' > sequence.txt; echo "7" >> sequence.txt'
    # tocall1 = 'cd ' + path + ' ; ' + tocall1
    os.system(tocall1)

    tocall2 = ' ./zagg'
    os.system(tocall2)

    res = np.loadtxt('zagg7.txt')

    os.chdir(path0)

    return res


def catGRANULE(seq, proteinName='PROTEIN'):
    path0 = os.getcwd()
    path = '/Users/mmonti/Documents/Scienza/Mutations_and_functions/gianAlgorithms/catGRANULE'

    path_seq = '/Users/mmonti/Documents/Scienza/Mutations_and_functions/gianAlgorithms/catGRANULE/sequences'
    os.chdir(path_seq)

    seq_file = open("sequence.txt", "w")
    seq_file.write('%s ' % (proteinName))

    seq_file.write('%s' % (seq))

    seq_file.close()

    os.chdir(path)

    tocall1 = 'sh run_catGRANULE.sh ' + ' > results.txt'

    # tocall1 = 'cd ' + path + ' ; ' + tocall1
    os.system(tocall1)
    dat = pd.read_csv('results.txt', sep=' ')
    scoreS = dat.values.tolist()
    score = float(scoreS[0][1])

    profile = np.array([[float(i[0]), float(i[1])] for i in scoreS[2:]])

    os.chdir(path0)

    return score, profile


def Zscore(a):
    return ((a - np.mean(a)) / np.std(a))


def DeltaZyggScore(seq_wt, mut_seq):
    zw = Ziggregator(seq_wt)[0]
    cgw = catGRANULE(seq_wt)[0]
    cfw = camFOLD(seq_wt)[0]
    proZw = np.array(Ziggregator_profile(seq_wt)[:, 1])
    proGw = np.array(catGRANULE(seq_wt)[1][:, 1])
    proFw = np.array(camFOLD_profile(seq_wt)[1]) + 0.000001

    zm = Ziggregator(mut_seq)[0]
    cgm = catGRANULE(mut_seq)[0]
    cfm = camFOLD(mut_seq)[0]
    proZm = np.array(Ziggregator_profile(mut_seq)[:, 1])
    proGm = np.array(catGRANULE(mut_seq)[1][:, 1])
    proFm = np.array(camFOLD_profile(mut_seq)[1]) + 0.000001

    # simple difference between WT and mutSeq
    score1 = [zw - zm, cgw - cgm, cfw - cfm]

    # percentage difference between WT and mutSeq
    score2 = [(zw - zm) / zw, (cgw - cgm) / cgm, (cfw - cfm) / cfm]

    # from the profiles Mean((wt-mut)/wt)
    score3 = [np.mean((proZw - proZm) / proZw), np.mean((proGw - proGm) / proGw), np.mean((proFw - proFm) / proFw)]

    # from the profiles Mean((wt-mut)/(wt+mut))
    score4 = [np.mean((proZw - proZm) / (proZw + proZm)), np.mean((proGw - proGm) / (proGw + proGm)),
              np.mean((proFw - proFm) / (proFw + proFm))]

    return score1, score2, score3, score4


def DeltaZyggScore2(seq_wt, mut_seq):
    proZw = np.array(Ziggregator_profile(seq_wt)[:, 1])
    proGw = np.array(catGRANULE(seq_wt)[1][:, 1])
    proFw = np.array(camFOLD_profile(seq_wt)[1]) + 0.000001

    zw = np.mean(proZw)
    cgw = np.mean(proGw)
    cfw = np.mean(proFw)

    proZm = np.array(Ziggregator_profile(mut_seq)[:, 1])
    proGm = np.array(catGRANULE(mut_seq)[1][:, 1])
    proFm = np.array(camFOLD_profile(mut_seq)[1]) + 0.000001

    zm = np.mean(proZm)
    cgm = np.mean(proGm)
    cfm = np.mean(proFm)

    # simple difference between WT and mutSeq
    score1 = [zw - zm, cgw - cgm, cfw - cfm]

    # percentage difference between WT and mutSeq
    score2 = [(zw - zm) / zw, (cgw - cgm) / cgm, (cfw - cfm) / cfm]

    # from the profiles Mean((wt-mut)/wt)
    score3 = [np.mean((proZw - proZm) / proZw), np.mean((proGw - proGm) / proGw), np.mean((proFw - proFm) / proFw)]

    # from the profiles Mean((wt-mut)/(wt+mut))
    score4 = [np.mean((proZw - proZm) / (proZw + proZm)), np.mean((proGw - proGm) / (proGw + proGm)),
              np.mean((proFw - proFm) / (proFw + proFm))]

    return score1, score2, score3, score4


def printListOnString(itemlist, separator="\t"):
    string = separator.join(str(item) for item in itemlist)
    return string


def mutationPiece(seq_wt, mutPiece, start=0):
    seq_wt = str(seq_wt)
    newSeq = seq_wt[:start] + mutPiece + seq_wt[start + len(mutPiece):]

    return newSeq


def translateSeqtoNumbers(seq):
    letters = seq

    numbers = []
    for letter in letters:
        number = ord(letter) - 64
        numbers.append(number)

    return numbers


def zygFoldGran_correlation(seq_wt):
    profileWT_c2 = catGRANULE(seq_wt)[1][:, 1]

    profileWT_Zigg = Ziggregator_profile(seq_wt)[:, 1]

    profile_camFOLD = camFOLD_profile(seq_wt)[1]

    Zyg = smooth2(profileWT_Zigg, box_pts=51)
    catG = smooth2(profileWT_c2, box_pts=1)
    camF = smooth2(profile_camFOLD, box_pts=45)

    cor1, pv1 = pearsonr(Zyg, catG)

    cor2, pv2 = pearsonr(camF, catG)

    cor3, pv3 = pearsonr(Zyg, camF)

    # print(len(profileWT_Zigg),len(profileWT_c2),len(profile_camFOLD))
    return cor1, cor2, cor3


def zygFoldGran_correlation_2(seq_wt):
    profileWT_c2 = catGRANULE(seq_wt)[1][:, 1]

    profileWT_Zigg = Ziggregator_profile(seq_wt)[:, 1]

    profile_camFOLD = camFOLD_profile(seq_wt)[1]

    Zyg = smooth2(profileWT_Zigg, box_pts=51)
    catG = smooth2(profileWT_c2, box_pts=1)
    camF = smooth2(profile_camFOLD, box_pts=45)

    cor1, pv1 = pearsonr(Zyg, catG)

    cor2, pv2 = pearsonr(camF, catG)

    cor3, pv3 = pearsonr(Zyg, camF)

    # print(len(profileWT_Zigg),len(profileWT_c2),len(profile_camFOLD))
    return cor1, cor2, cor3, Zyg, catG, camF


def fromUniprot(cID):
    # better to use this code and loop from script because uniprot could respond slowly and overcharger the server
    baseUrl = "http://www.uniprot.org/uniprot/"
    currentUrl = baseUrl + cID + ".fasta"
    response = r.post(currentUrl)
    cData = ''.join(response.text)

    Seq = StringIO(cData)
    pSeq = list(SeqIO.parse(Seq, 'fasta'))

    return pSeq


def fromUniprot_list(ListID):
    baseUrl = "http://www.uniprot.org/uniprot/"
    pSeq = []
    for cID in ListID:
        currentUrl = baseUrl + cID + ".fasta"
        response = r.post(currentUrl)
        cData = ''.join(response.text)

        Seq = StringIO(cData)
        pSeq.append(list(SeqIO.parse(Seq, 'fasta'))[0])

    return pSeq


def common_elements(list1, list2, unique=False):
    result = []
    pos = []
    count = 0
    for element in list1:
        if element in list2:
            result.append(element)
            pos.append(count)
        count = count + 1
    if unique:
        return np.unique(result), pos

    return result, pos


def shortScale(a):
    b = (np.zeros(len(a) - 6) + 1)
    d = np.concatenate((np.array([7, 7. / 3., 7. / 5.]), b, np.array([7. / 5., 7. / 3., 7.])), axis=None)
    print(len(a), len(b))
    return a * d




def compute_contact_map(pdb_file, cutoff=6.):
    parser = PDBParser()
    structure = parser.get_structure('pdb', pdb_file)
    ppb = PPBuilder()
    model = structure[0]

    residues = list(model.get_residues())
    n_residues = len(residues)
    contact_map = np.zeros((n_residues, n_residues))

    for i in range(n_residues):
        for j in range(i + 1, n_residues):
            residue_i = residues[i]
            residue_j = residues[j]
            if residue_i.is_disordered() or residue_j.is_disordered():
                continue
            if residue_i.has_id('CA') and residue_j.has_id('CA'):
                ca_distance = residue_i['CA'] - residue_j['CA']
                if ca_distance < cutoff:
                    contact_map[i][j] = 1
                    contact_map[j][i] = 1

    return contact_map


def compute_protein_size(pdb_file):
    # it is like length of the sequence so it could be eliminate
    parser = PDBParser()
    structure = parser.get_structure('pdb', pdb_file)
    model = structure[0]
    residues = list(model.get_residues())
    n_residues = len(residues)
    return n_residues


def computeSize_GyrationRadius(pdb_file):
    # Load PDB file
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)

    # Compute center of mass
    masses = []
    positions = []
    for atom in structure.get_atoms():
        element = atom.element.strip()
        if element in ['C', 'N', 'O', 'S']:
            mass = atom.mass
            position = atom.get_coord()
            masses.append(mass)
            positions.append(position)
    masses = np.array(masses)
    positions = np.array(positions)
    center_of_mass = np.sum(positions * masses[:, np.newaxis], axis=0) / np.sum(masses)

    # Compute Rg
    residues = Selection.unfold_entities(structure, 'R')
    rg_sum = 0
    for residue in residues:
        position = residue['CA'].get_coord()
        rg_sum += np.sum((position - center_of_mass) ** 2)
    Rg = np.sqrt(rg_sum / len(residues))

    # print('Radius of gyration:', Rg)

    # Compute average distance among residues
    residues = Selection.unfold_entities(structure, 'R')
    distances = []
    for i, residue1 in enumerate(residues):
        for residue2 in residues[i + 1:]:
            distance = residue1['CA'] - residue2['CA']
            distances.append(distance)
    average_distance = sum(distances) / len(distances)

    # print('Average distance among residues:', average_distance)

    return Rg, average_distance


def compute_disorder_propensity(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('pdb', pdb_file)
    model = structure[0]
    residues = list(model.get_residues())
    disorder_scores = []
    for residue in residues:
        if residue.is_disordered():
            disorder_scores.append(1.0)
        else:
            disorder_scores.append(0.0)
    return (disorder_scores)


def getSequecePDB(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('pdb', pdb_file)
    ppb = PPBuilder()
    for pp in ppb.build_peptides(structure):
        (pp.get_sequence())
    seq = str(pp.get_sequence().__str__())
    return seq


def calculate_secondary_structure(pdb_file, plot=False):
    # step 1 define the sequence and the pdbparse object

    parser = PDBParser()
    structure = parser.get_structure('pdb_file', pdb_file)
    seqtr = getSequecePDB(pdb_file)

    # Step 2: Use the PPBuilder to build polypeptides
    ppb = PPBuilder()
    polypeptides = ppb.build_peptides(structure)

    # Step 3: Initialize the DSSP object and calculate secondary structure
    model = structure[0]
    dssp = DSSP(model, pdb_file, dssp='mkdssp')

    # step 4 iterate over DSSP object to compute the percentage of elements pof secondary structure

    a_key = list(dssp.keys())
    aminStruct = []
    for i in a_key:
        am = dssp[i]
        aminStruct.append(am[2])

    letter_counts = Counter(aminStruct)

    sectStructLetters = '-HBEGITS'
    letter_counts = [letter_counts[i] for i in sectStructLetters]

    if plot:
        plt.bar(list(sectStructLetters), height=np.array(letter_counts) / len(seqtr))

    # Step 6: Return the results as a array
    return list(np.array(letter_counts) / len(seqtr))


def compute_exposed_residues(pdb_file, threshold=0.5):
    parser = PDBParser()
    structure = parser.get_structure('pdb_file', pdb_file)

    model = structure[0]

    dssp = DSSP(model, pdb_file, dssp='mkdssp')

    a_key = list(dssp.keys())
    asa = []
    for i in a_key:
        am = dssp[i]
        asa.append(am[3])

    expRes = [1 if i > threshold else 0 for i in asa]

    return asa, expRes


def exposed_sequence(pdb_file, threshold=0.5):
    seqtr = getSequecePDB(pdb_file)

    asa, resEx = compute_exposed_residues(pdb_file, threshold)

    seqExp = [seqtr[i] for i in range(len(seqtr)) if resEx[i] == 1]
    seqInt = [seqtr[i] for i in range(len(seqtr)) if resEx[i] == 0]
    indExp = [i for i in range(len(seqtr)) if resEx[i] == 1]
    indInt = [i for i in range(len(seqtr)) if resEx[i] == 0]

    return ''.join(seqExp), ''.join(seqInt), indExp, indInt


def exposed_residues_properties(pdb_file, threshold=0.5):
    parser = PDBParser()
    structure = parser.get_structure('pdb_file', pdb_file)
    seqtr = getSequecePDB(pdb_file)

    model = structure[0]
    dssp = DSSP(model, pdb_file, dssp='mkdssp')

    a_key = list(dssp.keys())
    asa = []
    aminStruct = []
    for i in a_key:
        am = dssp[i]
        asa.append(am[3])
        aminStruct.append(am[2])

    expRes = [1 if i > threshold else 0 for i in asa]
    aminStruct_exp = [aminStruct[i] for i in range(len(seqtr)) if expRes[i] == 1]
    letter_counts = Counter(aminStruct_exp)
    sectStructLetters = '-HBEGITS'
    letter_counts = [letter_counts[i] for i in sectStructLetters]

    secStruct = np.array(letter_counts) / len(seqtr)

    return list(secStruct);


def histoletters(i, norm=True):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    letter_counts = Counter(i)
    letter_counts = [letter_counts[i] for i in amino_acids]

    if norm:
        compo = np.array(letter_counts) / len(i)
    else:
        compo = np.array(letter_counts)

    return compo


def sliding_window(string, window_size):
    # Convert string to NumPy array of characters
    char_array = np.array(list(string))
    shape = (len(string) - window_size + 1, window_size)
    strides = (char_array.strides[0], char_array.strides[0])
    windows = np.lib.stride_tricks.as_strided(char_array, shape=shape, strides=strides)
    windows = np.apply_along_axis(''.join, 1, windows)
    return windows


# apply the sliding window and compute the composition
def composition_slidingWindow(seqtr, window_sz=10):
    win = sliding_window(seqtr, window_sz)
    hiss = [histoletters(i) for i in win]

    return hiss


def slidingWindow_exposedMean(pdb_file, ws=10):
    # very similiar to :
    # expo2 ,pp=compute_exposed_residues(pdb_file)
    # expo2 = smooth(expo2,10)
    # but this takes into account exactly the position

    asa, res = compute_exposed_residues(pdb_file)
    seqtr = getSequecePDB(pdb_file)

    win = sliding_window(seqtr, ws)

    expo = [np.array(asa)[seqtr.find(i):seqtr.find(i) + ws].mean() for i in win]

    return expo


# for i in range(len(df)):
def compositionDIstance(seqtr, df, i, plot=False):
    comp = df['composition'][i]
    leng = df['length'][i]

    if len(seqtr) >= leng:

        hiss = composition_slidingWindow(seqtr, window_sz=leng)
        prof = [np.linalg.norm(i - comp) for i in hiss]
        mini = np.min(prof)
        medi = np.mean(prof)

    else:
        prof = np.empty((leng,))
        prof[:] = np.nan
        prof = list(prof)
        # print(prof)
        mini = np.nan
        medi = np.nan
    if plot:
        plt.plot(prof)
    return prof, mini, medi


# compute all the distance between the seqence sliding window and the composition vocabulary
def computeAllDistance_composition(seqtr, df):
    # seqtr = getSequecePDB(pdb_file)
    allcompo = []
    allmini = []
    allmedi = []

    for i in range(len(df)):
        a, b, c = compositionDIstance(seqtr, df, i)

        allcompo.append(a)
        allmini.append(b)
        allmedi.append(c)

    return allcompo, allmini, allmedi


# you can check where are the minima and are relatrive to which rbd in the database

def minimaPosition_bestFragmet(allmini, printout=False):
    minpos = np.where(allmini == np.min(allmini))[0]
    rbdfrag = df['peptide_chain'].values[minpos]
    if printout:
        print('minima positions in the df:', minpos)

    count = 0
    minrel = []
    for i in rbdfrag:
        if printout:
            print('reference fragment: %s, relative distance:' % i, np.array(allmini)[minpos[count]])
        minrel.append(np.array(allmini)[minpos[count]])
        count = count + 1

    if printout == False:
        return minpos, rbdfrag, minrel


def compute_charge(sequence):
    aa_charges = {
        'C': -0.04, 'D': -1.0, 'E': -1.0, 'H': 0.0, 'K': 1.0,
        'R': 1.0, 'Y': 0.0, 'A': 0.0, 'F': 0.0, 'G': 0.0,
        'I': 0.0, 'L': 0.0, 'M': 0.0, 'N': 0.0, 'P': 0.0,
        'Q': 0.0, 'S': 0.0, 'T': 0.0, 'V': 0.0, 'W': 0.0,
    }

    charge = sum([aa_charges[aa] for aa in sequence])
    return charge


## I ADDED THIS numba DECORATOR
# @cuda.jit
def pdbFeatures(pdb_file, df, thrASA=.5, cutoffCM=6):
    # exctract seuqence and internal and exposed part
    seqtr = getSequecePDB(pdb_file)
    seqExp, seqInt, indExp, indInt = exposed_sequence(pdb_file, threshold=thrASA)

    # compute the number of contacts, dimension, radius of giration, length, charge and
    # secondary structure features of external and full protein

    n_contaxct = np.sum(compute_contact_map(pdb_file, cutoff=cutoffCM))
    RG_protein, rmsd_prot = computeSize_GyrationRadius(pdb_file)
    plen = len(seqtr)
    secStructFull = calculate_secondary_structure(pdb_file, plot=False)
    secStruct_ext = exposed_residues_properties(pdb_file)
    fullCharge = compute_charge(seqtr)
    extCharge = compute_charge(seqExp)

    asa, rext = compute_exposed_residues(pdb_file, threshold=thrASA)

    # compute for a given sequence the RBD similarity of internal external and full sequence

    allcompo_ex, allmini_ex, allmedi_ex = computeAllDistance_composition(seqExp, df)
    allcompo_int, allmini_int, allmedi_int = computeAllDistance_composition(seqInt, df)
    allcompo, allmini, allmedi = computeAllDistance_composition(seqtr, df)

    results = [n_contaxct, RG_protein, rmsd_prot, plen, secStructFull, secStruct_ext, fullCharge, extCharge, asa,
               allmini_ex, allmini_int, allmini]

    return results


def get_alphafold_features(af_code_dir, af_model_dir):
    # In[6]:
    sn = 'RBDpep'
    # sn='CandidateRBDpep'
    # sn='Input'
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    aminoList = [i for i in amino_acids]

    rnadb = pd.read_excel('./src/AlphaFold/1-s2.0-S1097276516302878-mmc2.xlsx', sheet_name=sn)

    rbd = rnadb['MS-identified tryptic peptide (N-link)'].values
    uniid = rnadb['ProtID'].values
    complist = []
    lenlist = []
    for i in rbd[:]:
        letter_counts = Counter(i)
        letter_counts = [letter_counts[i] for i in amino_acids]
        # compo=np.array(list(letter_counts.values()))/len(i)
        compo = np.array(letter_counts) / len(i)

        # df = pd.DataFrame(compo,index=aminoList)
        # df.plot(kind='bar')
        complist.append(compo)
        lenlist.append(len(i))

    df = pd.DataFrame(list(zip(list(uniid), complist, lenlist, list(rbd))),
                      columns=['uniprotID', 'composition', 'length', 'peptide_chain'])


    t = .5  # percentage of available surface
    c = 6  # amstrong distance, if two residues are closer they are considered in contact

    # Read the names of the pdb files in the Alphafold folder
    pdbFileNames = glob(af_model_dir + '*pdb')
    # print(os.listdir('./alphafold_human_pdb/')[0])
    pdbFileNames = [i for i in pdbFileNames if 'pdb' in i]
    pdbIDs = [(os.path.basename(i)).split('.')[0] for i in pdbFileNames]

    X = [(f, df, t, c) for f in pdbFileNames]

    results = []
    n_cores = 22

    if n_cores == None or n_cores > mp.cpu_count():
        pool = mp.Pool(mp.cpu_count())
    else:
        pool = mp.Pool(n_cores)  # If the user wants to use a specified number of cores

    results = pool.starmap_async(pdbFeatures, X)
    results = results.get()
    pool.close()

    # print(type(results))
    # print(results)

    resdf = pd.DataFrame(results)
    uniprot_series = pd.Series(data=pdbIDs, name='Uniprot_ID')
    resdf.columns = ['n_contact', 'RG_protein', 'rmsd_prot', 'plen', 'secStructFull', 'secStruct_ext', 'fullCharge',
                     'extCharge', 'asa', 'allmini_ex', 'allmini_int', 'allmini']
    # resdf['Uniprot_ID']=pdbIDs
    resdf = pd.concat([uniprot_series, resdf], axis=1)
    resdf.to_pickle(out_folder + "alphafold_features.pkl")
    return out_folder + "alphafold_features.pkl"
    #resdf.to_csv(out_folder + "alphafold_features.csv", index=False)
    #return out_folder + "alphafold_features.csv"
    
def get_alphafold_features_from_file(af_code_dir, pdb_files):
    # In[6]:
    sn = 'RBDpep'
    # sn='CandidateRBDpep'
    # sn='Input'
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    aminoList = [i for i in amino_acids]

    rnadb = pd.read_excel('./src/AlphaFold/1-s2.0-S1097276516302878-mmc2.xlsx', sheet_name=sn)

    rbd = rnadb['MS-identified tryptic peptide (N-link)'].values
    uniid = rnadb['ProtID'].values
    complist = []
    lenlist = []
    for i in rbd[:]:
        letter_counts = Counter(i)
        letter_counts = [letter_counts[i] for i in amino_acids]
        # compo=np.array(list(letter_counts.values()))/len(i)
        compo = np.array(letter_counts) / len(i)

        # df = pd.DataFrame(compo,index=aminoList)
        # df.plot(kind='bar')
        complist.append(compo)
        lenlist.append(len(i))

    df = pd.DataFrame(list(zip(list(uniid), complist, lenlist, list(rbd))),
                      columns=['uniprotID', 'composition', 'length', 'peptide_chain'])


    t = .5  # percentage of available surface
    c = 6  # amstrong distance, if two residues are closer they are considered in contact

    # Read the names of the pdb files in the Alphafold folder
    #pdbFileNames = [pdb_file]
    pdbFileNames = pdb_files
    # print(os.listdir('./alphafold_human_pdb/')[0])
    pdbFileNames = [i for i in pdbFileNames if 'pdb' in i]
    pdbIDs = [(os.path.basename(i)).split('.')[0] for i in pdbFileNames]

    X = [(f, df, t, c) for f in pdbFileNames]

    results = []
    n_cores = 22

    if n_cores == None or n_cores > mp.cpu_count():
        pool = mp.Pool(mp.cpu_count())
    else:
        pool = mp.Pool(n_cores)  # If the user wants to use a specified number of cores


    results = pool.starmap_async(pdbFeatures, X)
    results = results.get()
    pool.close()

    resdf = pd.DataFrame(results)
    uniprot_series = pd.Series(data=pdbIDs, name='Uniprot_ID')
    resdf.columns = ['n_contact', 'RG_protein', 'rmsd_prot', 'plen', 'secStructFull', 'secStruct_ext', 'fullCharge',
                     'extCharge', 'asa', 'allmini_ex', 'allmini_int', 'allmini']
    # resdf['Uniprot_ID']=pdbIDs
    resdf = pd.concat([uniprot_series, resdf], axis=1)
    return resdf



def add_plddt_to_af_features(models_dir, af_features_file, output_dir):
    data = []
    pdbs = glob(models_dir + '*.pdb')
    #awk_path = '/mnt/c/Users/Light/Desktop/Italy/awk_get_plddt_scores.awk'
    awk_path = '/data01/dimitris/scripts/awk_get_plddt_scores.awk'
    #expr = r'(.*_\d+_[A-Z]+)_unrelaxed.*'
    # uncomment below if you downloaded for alphafold directly
    expr = r'AF-(.*)-F\d+.*'

    # for each pdb in the directory...
    for pdb in pdbs:
        # set the awk command that calculates the plddt scores
        cmd = ['awk', '-f', awk_path, pdb]
        # get the basename of the model
        bn = os.path.basename(pdb)
        m = re.match(expr, bn)
        # get the protein name
        protein_name = m.group(1)
        # execute the awk command
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # get the result and split it in lines
        stdout_list = [line for line in res.stdout.split('\n') if line]
        # convert the list to a numpy array
        np_stdout = np.array(stdout_list)
        # convert the plddt scores to floats
        np_stdout_f = np_stdout.astype(np.float64)
        # append to the empty list
        data.append([protein_name, np.mean(np_stdout_f), np.std(np_stdout_f), np_stdout_f])

    df = pd.DataFrame(data=data, columns=['name', 'plddt_mean', 'plddt_std_dev', 'plddt'])
    plddt_output_file_path = os.path.join(output_dir, 'proteins_plddt.csv')
    df.to_csv(plddt_output_file_path, index=False)

    # parse the alphafold features to a dataframe
    af_features_df = pd.read_pickle(af_features_file)
    #af_features_df = pd.read_csv(af_features_file)
    # get the correct name
    af_features_df['Uniprot_ID'] = af_features_df['Uniprot_ID'].apply(lambda x: re.match(expr, x).group(1))


    mean_plddts = []
    std_plddts = []
    for i, r in af_features_df.iterrows():
        # m = re.match(r'AF-(.*)-F\d+.*', r['Uniprot_ID'])
        # name = m.group(1)
        name = r['Uniprot_ID']
        mean_plddt = df.loc[df['name'] == name]['plddt_mean'].values[0]
        std_plddt = df.loc[df['name'] == name]['plddt_std_dev'].values[0]
        # print(mean_plddt, std_plddt)
        mean_plddts.append(mean_plddt)
        std_plddts.append(std_plddt)

    af_features_df['mean_plddt'] = mean_plddts
    af_features_df['std_plddt'] = std_plddts

    output_features_file_path = os.path.join(output_dir, 'alphafold_features_plddt.csv')
    af_features_df.to_csv(output_features_file_path, index=False)
    af_features_df.to_pickle(os.path.join(output_dir, 'alphafold_features_plddt.pkl'))
    return af_features_df
    

def add_plddt_to_af_features_from_file(pdb_files, af_features_df):
    data = []
    #pdbs = [pdb_file]
    pdbs = pdb_files
    awk_path = './src/awk_get_plddt_scores.awk'
    expr = r'AF-(.*)-F\d+.*'

    # for each pdb in the directory...
    for pdb in pdbs:
        # set the awk command that calculates the plddt scores
        cmd = ['awk', '-f', awk_path, pdb]
        # get the basename of the model
        #bn = os.path.basename(pdb)
        #m = re.match(expr, bn)
        # get the protein name
        #protein_name = m.group(1)
        # get the protein_name
        protein_name = os.path.basename(pdb).split('.')[0]
        # execute the awk command
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # get the result and split it in lines
        stdout_list = [line for line in res.stdout.split('\n') if line]
        # convert the list to a numpy array
        np_stdout = np.array(stdout_list)
        # convert the plddt scores to floats
        np_stdout_f = np_stdout.astype(np.float64)
        # append to the empty list
        data.append([protein_name, np.mean(np_stdout_f), np.std(np_stdout_f), np_stdout_f])

    df = pd.DataFrame(data=data, columns=['name', 'plddt_mean', 'plddt_std_dev', 'plddt'])
    # get the correct name | UNCOMMENT THIS IF YOU HAVE AF MODELS
    #af_features_df['Uniprot_ID'] = af_features_df['Uniprot_ID'].apply(lambda x: re.match(expr, x).group(1))

    mean_plddts = []
    std_plddts = []
    for i, r in af_features_df.iterrows():
        # m = re.match(r'AF-(.*)-F\d+.*', r['Uniprot_ID'])
        # name = m.group(1)
        name = r['Uniprot_ID']
        mean_plddt = df.loc[df['name'] == name]['plddt_mean'].values[0]
        std_plddt = df.loc[df['name'] == name]['plddt_std_dev'].values[0]
        # print(mean_plddt, std_plddt)
        mean_plddts.append(mean_plddt)
        std_plddts.append(std_plddt)

    af_features_df['mean_plddt'] = mean_plddts
    af_features_df['std_plddt'] = std_plddts

    return af_features_df


def transformSecDict(seq,dicto):
    # takes in a sequence seq  and a dictionary 'char_to_number'  and build up the list of numbers
    seq=list(seq)

    # Transform the characters into consecutive numbers
    return  [dicto[char] for char in seq]


def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def trasformSeqFromJson(seq,file_name, smth = 7):

    # take a json file where there is a dictiopnary of aminoacid substitution

    # Open the file for reading
    with open(file_name, "r") as file:
        # Load the JSON data from the file and parse it into a dictionary
        loaded_dict = json.load(file)

    # Now, 'loaded_dict' contains the dictionary loaded from the JSON file
    #print(loaded_dict)

    return smooth(np.array(transformSecDict(seq,loaded_dict)),smth)


def compute_chemphysProfiles(seq,listofscalesJSON, smth=1, chargeNorm=False):

    profiles=[]
    for i in listofscalesJSON:
        if chargeNorm:


            if 'charge' in i:

                profiles.append(smooth((trasformSeqFromJson(seq,i, smth=1)+1)/2,smth))
            else:

                profiles.append(trasformSeqFromJson(seq,i, smth=smth))
        else:

            profiles.append(trasformSeqFromJson(seq,i, smth=smth))

    return np.array(profiles)
    #return np.mean(np.array(profiles), axis=1)


def compute_chemphysProperties(seq,listofscalesJSON, smth=1, chargeNorm=False):

    profiles=[]
    for i in listofscalesJSON:
        if chargeNorm:


            if 'charge' in i:

                profiles.append(smooth((trasformSeqFromJson(seq,i, smth=1)+1)/2,smth))
            else:

                profiles.append(trasformSeqFromJson(seq,i, smth=smth))
        else:

            profiles.append(trasformSeqFromJson(seq,i, smth=smth))

    #return np.array(profiles), np.mean(np.array(profiles),axis=1)

    return np.mean(np.array(profiles), axis=1)


def get_list_from_pdb_string(pdb_files_string):
    if ',' in pdb_files_string:
        pdbFileNames = pdb_files_string.split(',')
        pdbFileNames = [x.strip() for x in pdbFileNames]
    else:
        pdbFileNames = pdb_files_string.split()
        
    if len(pdbFileNames) == 0:
        print("invalid pdb files")
        
    if any(not x.endswith('pdb') for x in pdbFileNames):
        print('you should provide valid pdb files')
        
    return pdbFileNames


def get_seq_from_pdb(pdb_file):
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    seq = []

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("prot", pdb_file)
    for res in structure.get_residues():
        seq.append(d3to1[res.get_resname()])
    return ''.join(seq)

  
    
def get_sequences_from_pdbs(pdbFileNames):
    #pdbFileNames = glob(models_dir + '*pdb')

    #pdbFileNames = [pdb_file]
    pdbIDs = [(os.path.basename(i)).split('.')[0] for i in pdbFileNames]

    n_cores = 22
    sequences = []

    if n_cores == None or n_cores > mp.cpu_count():
        pool = mp.Pool(mp.cpu_count())
    else:
        pool = mp.Pool(n_cores)  # If the user wants to use a specified number of cores

    sequences = pool.map_async(get_seq_from_pdb, pdbFileNames)
    sequences = sequences.get()
    pool.close()

    return sequences, pdbIDs

def get_sequences_from_aa_text(aa_text):
    names = []
    sequences = []
    names.append('Seq1')
    sequences.append(aa_text.replace('\n', ''))

    return sequences, names

def get_sequences_from_fasta(fasta_text):
    names = []
    sequences = []

    # Convert the input text to a file-like object
    text_file = StringIO(fasta_text)

    # Parse the sequences using SeqIO
    for record in SeqIO.parse(text_file, "fasta"):
        names.append(record.id)
        sequences.append(str(record.seq))

    return sequences, names


def get_sequences_from_file(fasta):
    names = []
    sequences = []

    # Parse the sequences using SeqIO
    for record in SeqIO.parse(fasta, "fasta"):
        names.append(record.id)
        sequences.append(str(record.seq))

    return sequences, names


def get_physical_chemical_properties(sequences, ids, scales_dir):
    sys.path.append(scales_dir)
    names = glob(f'{scales_dir}/*')
    names1 = [x.replace(f'{scales_dir}/', '') for x in names]
    names1 = [x.replace('.json', '') for x in names1]

    n_cores = 22

    X = [(f, names, 7, False) for f in sequences]

    results = []

    if n_cores == None or n_cores > mp.cpu_count():
        pool = mp.Pool(mp.cpu_count())
    else:
        pool = mp.Pool(n_cores)  # If the user wants to use a specified number of cores


    pl = pool.starmap_async(compute_chemphysProperties, X)
    results = pl.get()
    pool.close()

    resdf = pd.DataFrame(data=results, columns=names1, index=ids)
    return resdf


def get_physical_chemical_profiles(sequences, scales_dir, classifiers_dir, correct_order_columns):
    sys.path.append(scales_dir)
    names = glob(f'{scales_dir}/*')
    names1 = [x.replace(f'{scales_dir}/', '') for x in names]
    names1 = [x.replace('.json', '') for x in names1]

    n_cores = 22

    X = [(f, names, 1, False) for f in sequences]

    results = []

    if n_cores == None or n_cores > mp.cpu_count():
        pool = mp.Pool(mp.cpu_count())
    else:
        pool = mp.Pool(n_cores)  # If the user wants to use a specified number of cores

    pl = pool.starmap_async(compute_chemphysProfiles, X)
    results = pl.get()
    pool.close()
    
    rf_classifier = joblib.load(classifiers_dir + "ONLY_PHYSCHEM/RandomForest/gridsearchCV_Object.pkl")
    profiles = []
    # for each protein...
    for i in range(len(results)):
        # ...set the protein matrix
        mat = results[i]
        # create a dataframe with indexes the name of the features
        mat_df = pd.DataFrame(mat, index=names1)
        # ... reindex the dataframe in order to get the correct order
        correct_order_df = mat_df.reindex(index=correct_order_columns[:82])
        # ... compute the profile
        profile = ComputeProfile_fromMatrix2(correct_order_df, rf_classifier)
        # ... append it in the profiles array
        profiles.append(np.array(profile))
   
    return profiles


def predict(df, classifiers_dir, only_pc=False):
    '''
    A functions to calculate the predictions of both alphafold and physical/chemical properties
    and return a dictionary of the predictions
    '''
    classifier = 'MLP'
    all_feats_folder = './src/TRAINED_MODELS/ALL_FEATURES/'
    if only_pc == True:
        all_feats_folder = './src/TRAINED_MODELS/ONLY_PHYSCHEM/'
        classifier = 'RandomForest'

    dic = {}

    f2 = all_feats_folder + classifier + '/'

    model = joblib.load(f2 + "gridsearchCV_Object.pkl")

    pred = model.best_estimator_.predict_proba(df)[:, 1]

    dic[classifier] = pred

    return dic


def smooth_matrix_x(matrix, window=10):
    # Create a 1D smoothing kernel (e.g., Gaussian)
    kernel_size = window
    kernel = np.ones(kernel_size) / kernel_size

    # Apply convolution in the x-direction
    # smoothed_matrix = np.apply_along_axis(lambda x: convolve(x, kernel, mode='same'), axis=1, arr=matrix)
    smoothed_matrix = np.apply_along_axis(lambda x: np.convolve(x, kernel, mode='same'), axis=1, arr=matrix)

    return smoothed_matrix


def ComputeProfile_fromMatrix2(X, pipe, window=1):
    # X is the matrixes of lenght of the sequ times the number of features
    # pipe is the pipeline that scores
    # matrix where at each line you have the profile of each chem phys properties ( not the mean but the vlaue per aminoacid)
    # pipe is the RandomForest classifier only with chem phhys

    rolli2 = smooth_matrix_x(X, window=window)
    # newprof2=pipe.bestestimator.predict_proba(rolli2.T)[:,1]
    newprof2 = pipe.best_estimator_.predict_proba(rolli2.T)[:, 1]
    return newprof2


def convert_to_list(string_value):
    try:
        return eval(string_value)
    except (ValueError, SyntaxError):
        return np.nan


# allmini_ex, allmini_int, allmini: save the three minimum values and their index in the array
def mySortFunc(l):
    larr = np.array(l)
    min_values = np.sort(larr)[:3]
    min_indices = np.argsort(larr)[:3]

    return min_values, min_indices;


def retMin(df, col):
    df = df.copy()
    mymins0 = []
    mymins1 = []
    mymins2 = []
    myindexes0 = []
    myindexes1 = []
    myindexes2 = []
    for i in df.index:
        mins, indexes = mySortFunc(df.loc[i, col])
        mymins0.append(mins[0])
        mymins1.append(mins[1])
        mymins2.append(mins[2])
        myindexes0.append(indexes[0])
        myindexes1.append(indexes[1])
        myindexes2.append(indexes[2])

    df[col + '_min_0'] = mymins0
    df[col + '_min_1'] = mymins1
    df[col + '_min_2'] = mymins2
    df[col + '_indexmin_0'] = myindexes0
    df[col + '_indexmin_1'] = myindexes1
    df[col + '_indexmin_2'] = myindexes2

    return df;


def ProcessAF_Features(df):
    column_dict = dict(zip(['n_contact',
                       'RG_protein',
                       'rmsd_prot',
                       'plen',
                       'fullCharge',
                       'extCharge',
                       'secStructFull_0',
                       'secStructFull_1',
                       'secStructFull_2',
                       'secStructFull_3',
                       'secStructFull_4',
                       'secStructFull_5',
                       'secStructFull_6',
                       'secStructFull_7',
                       'secStruct_ext_0',
                       'secStruct_ext_1',
                       'secStruct_ext_2',
                       'secStruct_ext_3',
                       'secStruct_ext_4',
                       'secStruct_ext_5',
                       'secStruct_ext_6',
                       'secStruct_ext_7',
                       'asa_mean',
                       'asa_std',
                       'allmini_min_0',
                       'allmini_min_1',
                       'allmini_min_2',
                       'allmini_indexmin_0',
                       'allmini_indexmin_1',
                       'allmini_indexmin_2',
                       'allmini_int_min_0',
                       'allmini_int_min_1',
                       'allmini_int_min_2',
                       'allmini_int_indexmin_0',
                       'allmini_int_indexmin_1',
                       'allmini_int_indexmin_2',
                       'allmini_ex_min_0',
                       'allmini_ex_min_1',
                       'allmini_ex_min_2',
                       'allmini_ex_indexmin_0',
                       'allmini_ex_indexmin_1',
                       'allmini_ex_indexmin_2',
                       'mean_plddt',
                       'std_plddt'],
                      ['n_contacts',
                       'RG_protein',
                       'rmsd_prot',
                       'Length',
                       'fullCharge',
                       'extCharge',
                       'Percentage_Coil_FullSeq',
                       'Percentage_AlphaHelix_FullSeq',
                       'Percentage_BetaBridge_FullSeq',
                       'Percentage_Strand_FullSeq',
                       'Percentage_3_10Helix_FullSeq',
                       'Percentage_PiHelix_FullSeq',
                       'Percentage_Turn_FullSeq',
                       'Percentage_Bend_FullSeq',
                       'Percentage_Coil_ExtSeq',
                       'Percentage_AlphaHelix_ExtSeq',
                       'Percentage_BetaBridge_ExtSeq',
                       'Percentage_Strand_ExtSeq',
                       'Percentage_3_10Helix_ExtSeq',
                       'Percentage_PiHelix_ExtSeq',
                       'Percentage_Turn_ExtSeq',
                       'Percentage_Bend_ExtSeq',
                       'asa_mean',
                       'asa_std',
                       'RBD_Full_min_0',
                       'RBD_Full_min_1',
                       'RBD_Full_min_2',
                       'RBD_Full_indexmin_0',
                       'RBD_Full_indexmin_1',
                       'RBD_Full_indexmin_2',
                       'RBD_int_min_0',
                       'RBD_int_min_1',
                       'RBD_int_min_2',
                       'RBD_int_indexmin_0',
                       'RBD_int_indexmin_1',
                       'RBD_int_indexmin_2',
                       'RBD_ext_min_0',
                       'RBD_ext_min_1',
                       'RBD_ext_min_2',
                       'RBD_ext_indexmin_0',
                       'RBD_ext_indexmin_1',
                       'RBD_ext_indexmin_2',
                       'average_plddt',
                       'stddev_plddt']))
    df['asa'] = df['asa'].apply(convert_to_list)
    df['allmini_ex'] = df['allmini_ex'].apply(convert_to_list)
    df['allmini_int'] = df['allmini_int'].apply(convert_to_list)
    df['allmini'] = df['allmini'].apply(convert_to_list)
    df['secStructFull'] = df['secStructFull'].apply(convert_to_list)
    df['secStruct_ext'] = df['secStruct_ext'].apply(convert_to_list)

    # Expand the column corresponding to secondary structure
    expanded_columns_secStructFull = df['secStructFull'].apply(pd.Series)
    expanded_columns_secStructFull.columns = ['secStructFull_' + str(i) for i in
                                              range(len(expanded_columns_secStructFull.columns))]
    expanded_columns_secStruct_ext = df['secStruct_ext'].apply(pd.Series)
    expanded_columns_secStruct_ext.columns = ['secStruct_ext_' + str(i) for i in
                                              range(len(expanded_columns_secStruct_ext.columns))]

    # Drop columns corresponding to secondary structure from dataset2
    df = df.drop('secStructFull', axis=1)
    df = df.drop('secStruct_ext', axis=1)

    # Assign the expanded columns back to the original DataFrame
    df[expanded_columns_secStructFull.columns] = expanded_columns_secStructFull
    df[expanded_columns_secStruct_ext.columns] = expanded_columns_secStruct_ext

    # asa: mean and variance
    asa_means = []
    asa_std = []

    for i in df.index:
        asa_means.append(np.mean(df.loc[i, 'asa']))
        asa_std.append(np.std(df.loc[i, 'asa']))

    df["asa_mean"] = asa_means
    df["asa_std"] = asa_std
    df = df.drop('asa', axis=1)

    df = retMin(df, 'allmini')
    df = retMin(df, 'allmini_int')
    df = retMin(df, 'allmini_ex')

    df = df.drop('allmini', axis=1)
    df = df.drop('allmini_int', axis=1)
    df = df.drop('allmini_ex', axis=1)

    df = df.rename(columns=column_dict)

    df['n_contacts_norm'] = df['n_contacts'] / df['Length']
    df['RG_protein_norm'] = df['RG_protein'] / df['Length']

    return df;


def ProcessAF_Features_from_file(df):
    column_dict = dict(zip(['n_contact',
                       'RG_protein',
                       'rmsd_prot',
                       'plen',
                       'fullCharge',
                       'extCharge',
                       'secStructFull_0',
                       'secStructFull_1',
                       'secStructFull_2',
                       'secStructFull_3',
                       'secStructFull_4',
                       'secStructFull_5',
                       'secStructFull_6',
                       'secStructFull_7',
                       'secStruct_ext_0',
                       'secStruct_ext_1',
                       'secStruct_ext_2',
                       'secStruct_ext_3',
                       'secStruct_ext_4',
                       'secStruct_ext_5',
                       'secStruct_ext_6',
                       'secStruct_ext_7',
                       'asa_mean',
                       'asa_std',
                       'allmini_min_0',
                       'allmini_min_1',
                       'allmini_min_2',
                       'allmini_indexmin_0',
                       'allmini_indexmin_1',
                       'allmini_indexmin_2',
                       'allmini_int_min_0',
                       'allmini_int_min_1',
                       'allmini_int_min_2',
                       'allmini_int_indexmin_0',
                       'allmini_int_indexmin_1',
                       'allmini_int_indexmin_2',
                       'allmini_ex_min_0',
                       'allmini_ex_min_1',
                       'allmini_ex_min_2',
                       'allmini_ex_indexmin_0',
                       'allmini_ex_indexmin_1',
                       'allmini_ex_indexmin_2',
                       'mean_plddt',
                       'std_plddt'],
                      ['n_contacts',
                       'RG_protein',
                       'rmsd_prot',
                       'Length',
                       'fullCharge',
                       'extCharge',
                       'Percentage_Coil_FullSeq',
                       'Percentage_AlphaHelix_FullSeq',
                       'Percentage_BetaBridge_FullSeq',
                       'Percentage_Strand_FullSeq',
                       'Percentage_3_10Helix_FullSeq',
                       'Percentage_PiHelix_FullSeq',
                       'Percentage_Turn_FullSeq',
                       'Percentage_Bend_FullSeq',
                       'Percentage_Coil_ExtSeq',
                       'Percentage_AlphaHelix_ExtSeq',
                       'Percentage_BetaBridge_ExtSeq',
                       'Percentage_Strand_ExtSeq',
                       'Percentage_3_10Helix_ExtSeq',
                       'Percentage_PiHelix_ExtSeq',
                       'Percentage_Turn_ExtSeq',
                       'Percentage_Bend_ExtSeq',
                       'asa_mean',
                       'asa_std',
                       'RBD_Full_min_0',
                       'RBD_Full_min_1',
                       'RBD_Full_min_2',
                       'RBD_Full_indexmin_0',
                       'RBD_Full_indexmin_1',
                       'RBD_Full_indexmin_2',
                       'RBD_int_min_0',
                       'RBD_int_min_1',
                       'RBD_int_min_2',
                       'RBD_int_indexmin_0',
                       'RBD_int_indexmin_1',
                       'RBD_int_indexmin_2',
                       'RBD_ext_min_0',
                       'RBD_ext_min_1',
                       'RBD_ext_min_2',
                       'RBD_ext_indexmin_0',
                       'RBD_ext_indexmin_1',
                       'RBD_ext_indexmin_2',
                       'average_plddt',
                       'stddev_plddt']))
    
    # Expand the column corresponding to secondary structure
    expanded_columns_secStructFull = df['secStructFull'].apply(pd.Series)
    expanded_columns_secStructFull.columns = ['secStructFull_' + str(i) for i in
                                              range(len(expanded_columns_secStructFull.columns))]
    expanded_columns_secStruct_ext = df['secStruct_ext'].apply(pd.Series)
    expanded_columns_secStruct_ext.columns = ['secStruct_ext_' + str(i) for i in
                                              range(len(expanded_columns_secStruct_ext.columns))]

    # Drop columns corresponding to secondary structure from dataset2
    df = df.drop('secStructFull', axis=1)
    df = df.drop('secStruct_ext', axis=1)

    # Assign the expanded columns back to the original DataFrame
    df[expanded_columns_secStructFull.columns] = expanded_columns_secStructFull
    df[expanded_columns_secStruct_ext.columns] = expanded_columns_secStruct_ext

    # asa: mean and variance
    asa_means = []
    asa_std = []

    for i in df.index:
        asa_means.append(np.mean(df.loc[i, 'asa']))
        asa_std.append(np.std(df.loc[i, 'asa']))

    df["asa_mean"] = asa_means
    df["asa_std"] = asa_std
    df = df.drop('asa', axis=1)

    df = retMin(df, 'allmini')
    df = retMin(df, 'allmini_int')
    df = retMin(df, 'allmini_ex')

    df = df.drop('allmini', axis=1)
    df = df.drop('allmini_int', axis=1)
    df = df.drop('allmini_ex', axis=1)

    df = df.rename(columns=column_dict)

    df['n_contacts_norm'] = df['n_contacts'] / df['Length']
    df['RG_protein_norm'] = df['RG_protein'] / df['Length']

    return df;


def untar(tar_file, models_dir):
    os.system(f'mkdir {models_dir}')
    cmd = f'tar -xvf {tar_file} -C {models_dir}'
    res = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(f'something went wrong in the untaring: {res.stderr}')
    else:
        print(f'untarring {os.path.basename(tar_file)} done')
    os.system(f'rm {models_dir}/*.cif.gz')
    
    
def return_af_protein_name(name):
    expr = r'AF-(.*)-F\d+.*'
    m = re.match(expr, name)
    # get the protein name
    protein_name = m.group(1)
    return protein_name
