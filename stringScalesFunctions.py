import numpy as np
from matplotlib import pyplot as plt
from Bio.Seq import Seq

from glob import glob

import math
from sklearn.decomposition import PCA


#from dtw import dtw
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
from scipy import interp
from sklearn.metrics import roc_auc_score

# Import some data to play with
from scipy.stats import pearsonr
import requests as r
from io import StringIO

from collections import Counter

def chunkstring(string, length):
    return [string[0+i:length+i] for i in range(0, len(string), length)]


def Stringentropy(string):
    "Calculates the Shannon entropy of a string"

    # get probability of chars in string
    prob = [ float(string.count(c)) / len(string) for c in dict.fromkeys(list(string)) ]

    # calculate the entropy
    entropy = - sum([ p * math.log(p) / math.log(2.0) for p in prob ])

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
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def smooth2(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='valid')
    return y_smooth




def remove_values_from_list(the_list, val):
    return [value for value in the_list if value != val]

def condensateMixtList(expT):
    
    return [np.concatenate( expT[i], axis=0 ) for i in range(np.size(expT))]


manhattan_distance = lambda x, y: np.abs(x - y)


def computeDeltaChemPhysPCA(ScaleNames, PCAcomp,WTseq, mutSeq):
    first =PCAcomp
    count = 0
    newScaleWT=[]

    newScaleMut=[]

    for fn in ScaleNames:
        newScaleWT.append(first[count]*np.array(transformSeq(fn, WTseq)))

        newScaleMut.append(first[count]*np.array(transformSeq(fn, mutSeq)))

        count = count+1
    amin_1 = np.array(newScaleWT)

        
    amin_2 = np.array(newScaleMut)
    
    return np.mean(amin_1-amin_2)


def computeMutationScorePCA(ScaleNames, PCAcomp,WTseq,mutation):

    numberMutations = np.shape(mutation)[0]
    mutScore = []
    bindingScore=[]
    for i in range(numberMutations):
        #print(i)

        index, wt, mut = mutation.iat[i,0],mutation.iat[i,1],mutation.iat[i,2]
        bindingScore.append(mutation.iat[i,3])
        mutSeq = seqMutation(WTseq, index, str(mut))

        deltaMut = computeDeltaChemPhysPCA(ScaleNames, PCAcomp,WTseq, mutSeq)
        mutScore.append(deltaMut)
        
    return bindingScore, mutScore
    
    
    


def transformSeq(filename, seq):

    fileT = open('seqTemporanea.txt', 'w')
    fileT.write(seq)
#     np.savetxt('seqTemporanea.txt',seq)
    fileT.close()


    tocall = 'bash ' + filename + ' seqTemporanea.txt -> TemporaneaChemPhys.txt' 
    os.system(tocall)

    #fileT2= open('TemporaneaChemPhys.txt','r')

    #ChemPhys_Seq = [i for i in fileT2.readlines()]

    ChemPhys_Seq = np.loadtxt('TemporaneaChemPhys.txt')
    os.system('rm seqTemporanea.txt')

    os.system('rm TemporaneaChemPhys.txt')

    return ChemPhys_Seq

def computeDeltaChemPhys(filename,WTseq, mutSeq):
    
    amin_1 = transformSeq(filename, WTseq)[:,1]

    amin_2 = transformSeq(filename, mutSeq)[:,1]

    
    return np.mean(amin_1-amin_2)
    
        
def seqMutation(seq, mutationIndex, mutatedAmino):


    mutSeq = seq[:mutationIndex] + mutatedAmino + seq[mutationIndex + 1:]
    return mutSeq



    
def computeMutationScore(filename,WTseq,mutation):
    fn = filename
    numberMutations = np.shape(mutation)[0]
    mutScore = []
    bindingScore=[]
    for i in range(numberMutations):
        #print(i)

        index, wt, mut = mutation.iat[i,0],mutation.iat[i,1],mutation.iat[i,2]
        bindingScore.append(mutation.iat[i,3])
        mutSeq = seqMutation(WTseq, index, str(mut))

        deltaMut = computeDeltaChemPhys(fn,WTseq, mutSeq)
        mutScore.append(deltaMut)
        
    return bindingScore, mutScore

def giveIndx(mutScore, tresh):
    return [(np.sign(i-tresh)+1)/2 for i in mutScore]

def diff_letters(a,b):
    return sum ( a[i] != b[i] for i in range(len(a)) )

def translate_from_dict(original_text,dictionary_of_translations):
    out = original_text
    for target in dictionary_of_translations:
        trans = str.maketrans(target,dictionary_of_translations[target])
        out = out.translate(trans)
    return out



def flipIndx(arr):
    return [1 if i==0 else 0 for i in arr]

aminoList = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

def RandomMix(seq,perc):
    seq=chunkstring(seq,1)
    seqNew=[]
    for i in seq:
        rn = np.random.random()
        if perc>rn:
            j=aminoList[np.random.randint(len(aminoList))]
        else:
            j=i
            
        seqNew.append(j)
        
    ns=''.join(seqNew)
    return ns
    
        

def multipleRandomSeq(seq,N,perc):
    
    seqList=[]
    for i in range(N):
        seqList.append(RandomMix(seq,perc))
        
    return seqList
        
def multipleRandomSeq_transform(fn,seq,N,perc):

    seqList=multipleRandomSeq(seq,N,perc)
    scaleList=[]
    for j in range(len(seqList)):
        scaledS = transformSeq(fn, seqList[j])
        scaleList.append(scaledS)
    return scaleList



def multipleRandomSeq_transformMeanScore(fn,seq,N,perc):

    seqList=multipleRandomSeq(seq,N,perc)
    scaleList=[]
    for j in range(len(seqList)):
        scaledS = transformSeq(fn, seqList[j])
        scaleList.append(np.mean(scaledS))
    return scaleList


def TopBotRND(aggScore,scoreList,fn,WtSeq_utile,perc, plot=False):
    
    sortAgg = np.sort(aggScore)
    indx_sort = np.argsort(aggScore)
    scorSort = np.array(scoreList)[indx_sort]

    topBotList =[10,200,500,2000,5000, int(np.size(sortAgg)/2)]
    
    topBotList =[10,30,50,70,100,150,200,300,400,500, int(np.size(sortAgg)/2)]


    roc_auc=[]

    #perc =0.05
    if plot:
        f,ax =plt.subplots(1,2, figsize=(30,30))


    for i in topBotList:

        top=sortAgg[:i]
        bot=multipleRandomSeq_transformMeanScore(fn,WtSeq_utile,i,perc)

        topScore = scorSort[:i]
        botScore = np.zeros(len(scorSort[-i:]))-1

        newAgg = np.concatenate([np.array(top),np.array(bot)])
        newScore = np.concatenate([np.array(topScore),np.array(botScore)])

        indx_agg = flipIndx(giveIndx(newScore,0))


        x,y,q= roc_curve(indx_agg, newAgg)
        auc_val = roc_auc_score(indx_agg, newAgg)
        roc_auc.append(auc_val)

        if plot:
            ax[0].plot(x, y, label='AUC = %0.2f  top-bottom = %.4f'%(auc_val, i))
            ax[0].legend(loc='best')

    if plot:
        ax[1].plot(topBotList,roc_auc)
        ax[1].set_xlabel('Top-Bottom size')
        ax[1].set_ylabel('AUC')
        plt.show()
        
    return roc_auc,topBotList


#indx on experimental
def TopBot(aggScore,scoreList,fn,WtSeq_utile, plot=False):
    '''
    indx onexperimental
    top and bottom on theory
    '''

    sortAgg = np.sort(aggScore)
    indx_sort = np.argsort(aggScore)
    scorSort = np.array(scoreList)[indx_sort]

    topBotList =[10,20,30,40,50,60,70,80,90,100,125,150,175,200,300,500,600,750,1000,1500,2000,3000,4500,5000,7500,10000,12000,15000, int(np.size(sortAgg)/2)]

    topBotList =[10,30,50,70,100,150,200,300,400,500, int(np.size(sortAgg)/2)]

    roc_auc=[]
    if plot:
        f,ax =plt.subplots(1,2, figsize=(30,30))


    for i in topBotList:

        top=sortAgg[:i]
        bot=sortAgg[-i:]

        topScore = scorSort[:i]
        botScore = scorSort[-i:]

        newAgg = np.concatenate([np.array(top),np.array(bot)])
        newScore = np.concatenate([np.array(topScore),np.array(botScore)])

        indx_agg = flipIndx(giveIndx(newScore,0))


        x,y,q= roc_curve(indx_agg, newAgg)
        auc_val = roc_auc_score(indx_agg, newAgg)
        roc_auc.append(auc_val)

        if plot:
            ax[0].plot(x, y, label='AUC = %0.2f  top-bottom = %.4f'%(auc_val, i))
            ax[0].legend(loc='best')

    if plot:
        ax[1].plot(topBotList,roc_auc,'.')

        ax[1].plot(topBotList,roc_auc,'-')

        ax[1].set_xlabel('Top-Bottom size')
        ax[1].set_ylabel('AUC')
        ax[1].set_xscale('log')



    return roc_auc,topBotList




#indx on experimental
def TopBot2(aggScore,scoreList,fn,WtSeq_utile, plot=False):
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


    topBotList =[10,20,30,40,50,60,70,80,90,100,125,150,175,200,300,500,600,750,1000,1500,2000,3000,4500,5000,7500,10000,12000,15000, int(np.size(sortAgg)/2)]

    #topBotList =[100,500,1000,5000,15000, int(np.size(sortAgg)/2)]

    roc_auc=[]
    if plot:
        f,ax =plt.subplots(1,2, figsize=(30,30))


    for i in topBotList:

        top=sortAgg[:i]
        bot=sortAgg[-i:]

        topScore = scorSort[:i]
        botScore = scorSort[-i:]

        newAgg = np.concatenate([np.array(top),np.array(bot)])
        newScore = np.concatenate([np.array(topScore),np.array(botScore)])

        indx_agg = flipIndx(giveIndx(newScore,0))


        x,y,q= roc_curve(indx_agg, newAgg)
        auc_val = roc_auc_score(indx_agg, newAgg)
        roc_auc.append(auc_val)

        if plot:
            ax[0].plot(x, y, label='AUC = %0.2f  top-bottom = %.4f'%(auc_val, i))
            ax[0].legend(loc='best')

    if plot:
        ax[1].plot(topBotList,roc_auc,'.')

        ax[1].plot(topBotList,roc_auc,'-')

        ax[1].set_xlabel('Top-Bottom size')
        ax[1].set_ylabel('AUC')
        ax[1].set_xscale('log')



    return roc_auc,topBotList


def dataVsRND(aggScore,sizeList,fn,WtSeq_utile,perc, plot=False):
    


    roc_auc=[]

    #perc =0.05
    if plot:
        f,ax =plt.subplots(1,2, figsize=(30,30))


    for i in sizeList:

        sortAgg=random.sample(aggScore,i)
        top=sortAgg
        bot=multipleRandomSeq_transformMeanScore(fn,WtSeq_utile,i,perc)

        topScore = np.zeros(len(sortAgg))+1
        botScore = np.zeros(len(sortAgg))

        newAgg = np.concatenate([np.array(top),np.array(bot)])
        indx_agg = np.concatenate([np.array(topScore),np.array(botScore)])

        #indx_agg = flipIndx(giveIndx(newScore,0))


        x,y,q= roc_curve(indx_agg, newAgg)
        auc_val = roc_auc_score(indx_agg, newAgg)
        roc_auc.append(auc_val)

        if plot:
            ax[0].plot(x, y, label='AUC = %0.2f  Pool size = %.4f'%(auc_val, i))
            ax[0].legend(loc='best')

    if plot:
        ax[1].plot(sizeList,roc_auc)
        ax[1].set_xlabel('Pool size')
        ax[1].set_ylabel('AUC')
        plt.show()
        
    return roc_auc,sizeList




def RandomMix_num(seq,num):
    seq=chunkstring(seq,1)
    seqNew=[]
    for i in range(num):
        rn = np.random.randint(len(seq))
        
        seq[rn]=aminoList[np.random.randint(len(aminoList))]
            
        
        
    ns=''.join(seq)
    return ns
    
    

def RandomMix_num_2(seq,num):
    seq=chunkstring(seq,1)
    seqNew=[]
    newAmList=[]
    oldAmList=[]
    positionMutation=[]
    for i in range(num):
        rn = np.random.randint(len(seq))
        
        newAm=aminoList[np.random.randint(len(aminoList))]
        oldAmList.append(seq[rn])
        newAmList.append(newAm)
        positionMutation.append(rn)
        
        seq[rn]=newAm
            
        
        
    ns=''.join(seq)
    return ns,oldAmList,newAmList,positionMutation
    
    
        

def multipleRandomSeq_num(seq,N,num):
    
    seqList=[]
    for i in range(N):
        seqList.append(RandomMix_num(seq,num))
        
    return seqList
        
def multipleRandomSeq_transform_num(fn,seq,N,num):

    seqList=multipleRandomSeq_num(seq,N,num)
    scaleList=[]
    for j in range(len(seqList)):
        scaledS = transformSeq(fn, seqList[j])
        scaleList.append(scaledS)
    return scaleList



def multipleRandomSeq_transformMeanScore_num(fn,seq,N,num):

    seqList=multipleRandomSeq_num(seq,N,num)
    scaleList=[]
    for j in range(len(seqList)):
        scaledS = transformSeq(fn, seqList[j])
        scaleList.append(np.mean(scaledS))
    return scaleList


def dataVsRND_num(aggScore,mutList,indxList,fn,WtSeq_utile, plot=False):
    


    roc_auc=[]

    #perc =0.05
    if plot:
        f,ax =plt.subplots(1,2, figsize=(30,30))


    for i in range(len(mutList)):
        print(len(indxList[i]))
        if(len(indxList[i]))>0:
            sortAgg=aggScore[indxList[i]]
            top=sortAgg
            bot=multipleRandomSeq_transformMeanScore_num(fn,WtSeq_utile,len(sortAgg),i)

            topScore = np.zeros(len(sortAgg))+1
            botScore = np.zeros(len(sortAgg))

            newAgg = np.concatenate([np.array(top),np.array(bot)])
            indx_agg = np.concatenate([np.array(topScore),np.array(botScore)])

            #indx_agg = flipIndx(giveIndx(newScore,0))


            x,y,q= roc_curve(indx_agg, newAgg)
            auc_val = roc_auc_score(indx_agg, newAgg)
            roc_auc.append(auc_val)

        if plot:
            ax[0].plot(x, y, label='AUC = %0.2f  N mutations = %.4f'%(auc_val, i))
            ax[0].legend(loc='best')

    if plot:
        ax[1].plot(range(len(mutList)),roc_auc)
        ax[1].set_xlabel('Pool size')
        ax[1].set_ylabel('AUC')
        plt.show()
        
    return roc_auc,range(len(mutList))






def dataVsRND_num_poolsize(aggScore,mutList,indxList,poolsize,fn,WtSeq_utile, plot=False):
    


    roc_auc=[]

    #perc =0.05
    if plot:
        f,ax =plt.subplots(1,2, figsize=(30,30))
    True_IL=[]
    
    for i in range(len(mutList)):
        if i%100==0:
            print(i/len(mutList))
        if(len(indxList[i]))>0:
            #print('NumberMutations:  ',i)


            sortAgg=aggScore[indxList[i]]
            if len(sortAgg)>poolsize:
                sortAgg=random.sample(list(sortAgg),poolsize)
            
            top=sortAgg
            bot=multipleRandomSeq_transformMeanScore_num(fn,WtSeq_utile,len(sortAgg),i)

            topScore = np.zeros(len(sortAgg))+1
            botScore = np.zeros(len(sortAgg))

            newAgg = np.concatenate([np.array(top),np.array(bot)])
            indx_agg = np.concatenate([np.array(topScore),np.array(botScore)])

            #indx_agg = flipIndx(giveIndx(newScore,0))


            x,y,q= roc_curve(indx_agg, newAgg)
            auc_val = roc_auc_score(indx_agg, newAgg)
            roc_auc.append(auc_val)
            True_IL.append(i)
        if plot:
            ax[0].plot(x, y, label='AUC = %0.2f  N mutations = %.d'%(auc_val, i))
            #ax[0].legend(loc='best')

    if plot:
        ax[1].plot(True_IL,roc_auc)
        ax[1].set_xlabel('Number of mutations')
        ax[1].set_ylabel('AUC')
        f.suptitle('Pool-Size %d'%poolsize)
        plt.show()
        
    return roc_auc,True_IL






def dataVsRND_num_poolsize0(aggScore,mutList,indxList,poolsize,fn,WtSeq_utile, plot=False):
    


    roc_auc=[]

    #perc =0.05
    if plot:
        f,ax =plt.subplots(1,2, figsize=(30,30))
    True_IL=[]
    
    for i in range(len(mutList)):
        if i%100==0:
            print(i/len(mutList))

        #print('NumberMutations:  ',i)


        sortAgg=aggScore[indxList[i]]
        if len(sortAgg)>poolsize:
            sortAgg=random.sample(list(sortAgg),poolsize)

        top=sortAgg
        bot=multipleRandomSeq_transformMeanScore_num(fn,WtSeq_utile,len(sortAgg),i)

        topScore = np.zeros(len(sortAgg))+1
        botScore = np.zeros(len(sortAgg))

        newAgg = np.concatenate([np.array(top),np.array(bot)])
        indx_agg = np.concatenate([np.array(topScore),np.array(botScore)])

        #indx_agg = flipIndx(giveIndx(newScore,0))


        x,y,q= roc_curve(indx_agg, newAgg)
        auc_val = roc_auc_score(indx_agg, newAgg)
        roc_auc.append(auc_val)
        True_IL.append(i)
        if plot:
            ax[0].plot(x, y, label='AUC = %0.2f  N mutations = %.d'%(auc_val, i))
            #ax[0].legend(loc='best')

    if plot:
        ax[1].plot(True_IL,roc_auc)
        ax[1].set_xlabel('Number of mutations')
        ax[1].set_ylabel('AUC')
        f.suptitle('Pool-Size %d'%poolsize)
        plt.show()
        
    return roc_auc,True_IL






def mutationDistribution(filename,WtSeq_utile):
    f = open(filename,'r')
    lines = f.readlines()





    difflet=[]

    for i in lines:
        mutSq=list(i[:-1])

        if len(mutSq)==len(WtSeq_utile):


            nl=diff_letters(WtSeq_utile,mutSq)
            #if nl<15:
             #   difflet.append(nl)

            difflet.append(nl)

    un = np.unique(difflet)

    mutList=[]
    indxList=[]
    for i in range(np.max(un)):

        mutList.append([])
        indxList.append([])


    count=0
    for i in lines:
        mutSq=list(i[:-1])

        if len(mutSq)==len(WtSeq_utile):
            nl=diff_letters(WtSeq_utile,mutSq)


            mutList[nl-1].append(mutSq)
            indxList[nl-1].append(count)

            count=count+1


    mutSize=[]
    for i in mutList:
        mutSize.append(len(i))
        
    return mutSize/np.sum(mutSize),mutList,indxList



#compute Right AggScore from sequences with same lenght

def aggSameLenght(aggScore, WtSeq_utile):
    aggScore2=[]
    aggScore3=[]
    for i in range(len(aggScore)):
        if len(WtSeq_utile)-1<seqSize[i]<=len(WtSeq_utile):
            aggScore2.append(aggScore[i])

        else :
            aggScore3.append(aggScore[i])

    print(len(aggScore),len(aggScore2),len(aggScore3))
    
    return aggScore2, aggScore3



def cutTresh(distr,tresh):
    
    return [1 if i> tresh else 0 for i in distr] 


#copmute agg score from fn scale and seqList list of sequences
def computeAggScore(seqList, fn):
    aggList=[]

    aggScore=[]

    for j in range(np.size(seqList)):

        if j%1000==0:
            print(j/np.size(seqList))

        scaledS = transformSeq(fn, seqList[j])
        aggList.append(scaledS)
        aggScore.append(np.mean(scaledS))
    
    return aggScore, aggList


def sequenceList(filename):
    

    f = open(filename,'r')
    lines = f.readlines()


    #generate the sequence list

    seqList =[]
    SeqId=[]
    count = 0
    for i in range(len(lines)):
        if i%2==1:
            seqList.append("".join(lines[i][:-1]))
        if i%2==0:
            SeqId.append(lines[i])
        count = count +1

    f.close()
    print(count)

    return seqList, SeqId

def LenghtDistribution(seqList, plot= False, density = True):


    seqSize=[]
    for i in seqList:
        seqSize.append(len(i))

    hist,bins2 =np.histogram(seqSize, bins=50, density = True)
    
    
    if plot:
        plt.plot(bins2[:-1],hist,'-')

        plt.xlabel('Sequence Lenght')
        plt.ylabel('Count')
        plt.show()
    
    return hist, bins2


def ScoreSeq(fn,seq):
    Amin0=transformSeq(fn, seq)[:,1]

    return np.mean(Amin0)


def div_d(my_dict):

    sum_p = sum(my_dict.values())

    for i in my_dict:
        my_dict[i] = float(my_dict[i]/sum_p)

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


def linearRandom(aminoDist,namesAminoOrdered):
    z= np.random.random()
    pre = 0
    post = 0 
    for i in range(0,len(aminoDist)):
        post = post +aminoDist[i]
        
        if pre<z<=post:

            return namesAminoOrdered[i]
        pre = pre + aminoDist[i]
        



def RandomMix_num_weighted(seq,num):
    seq=chunkstring(seq,1)
    seqNew=[]
    for i in range(num):
        rn = np.random.randint(len(seq))
        
        seq[rn]=linearRandom(aminoDist,namesAminoOrdered)
            
        
        
    ns=''.join(seq)
    return ns


def camFOLD(seq):
    path0=os.getcwd()
    path = '/Users/mmonti/Documents/Scienza/Mutations_and_functions/gianAlgorithms/camZFOLD/'

    os.chdir(path)

    tocall1 = 'sh start.sh ' + seq + ' > results.txt'

    #tocall1 = 'cd ' + path + ' ; ' + tocall1
    os.system(tocall1)
    f = open('results.txt', 'r')

    r = f.readlines()

    if r[2][44]=='[' or r[2][43]=='[':

        r1 = float(r[2][38:43])

        r2 = float(r[3][38:43])
    else:

        r1 = float(r[2][38:44])

        r2 = float(r[3][38:44])
    
    f.close()
    os.chdir(path0)


    return r1,r2

def camFOLD_profile(seq):
    path0=os.getcwd()
    path = '/Users/mmonti/Documents/Scienza/Mutations_and_functions/gianAlgorithms/camZFOLD/'

    os.chdir(path)

    tocall1 = 'sh start.sh ' + seq + ' > results.txt'

    #tocall1 = 'cd ' + path + ' ; ' + tocall1
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


    return x_indx , x_score


def Zscore(a):
    return((a-np.mean(a))/np.std(a))

def printListOnString(itemlist, separator = "\t"):
    
    string=separator.join(str(item) for item in itemlist)
    return string



def mutationPiece(seq_wt, mutPiece, start=0):
    seq_wt=str(seq_wt)
    newSeq = seq_wt[:start]+mutPiece+seq_wt[start+len(mutPiece):]
    
    return newSeq

def translateSeqtoNumbers(seq):
    
    letters = seq

    numbers = []
    for letter in letters:
        number = ord(letter) - 64
        numbers.append(number)

    return numbers
    
def fromUniprot(cID):
    #better to use this code and loop from script because uniprot could respond slowly and overcharger the server
    baseUrl="http://www.uniprot.org/uniprot/"
    currentUrl=baseUrl+cID+".fasta"
    response = r.post(currentUrl)
    cData=''.join(response.text)

    Seq=StringIO(cData)
    pSeq=list(SeqIO.parse(Seq,'fasta'))
    
    return pSeq



def fromUniprot_list(ListID):
    
    baseUrl="http://www.uniprot.org/uniprot/"
    pSeq=[]
    for cID in ListID:
        currentUrl=baseUrl+cID+".fasta"
        response = r.post(currentUrl)
        cData=''.join(response.text)

        Seq=StringIO(cData)
        pSeq.append(list(SeqIO.parse(Seq,'fasta'))[0])
    
    return pSeq


def common_elements(list1, list2, unique = False):
    result = []
    pos=[]
    count = 0
    for element in list1:
        if element in list2:
            result.append(element)
            pos.append(count)
        count=count+1
    if unique:
        return np.unique(result), pos
    
    return result,pos



def shortScale(a):

    b =(np.zeros(len(a)-6)+1)
    d=np.concatenate((np.array([7,7./3.,7./5.]), b,np.array([7./5.,7./3.,7.])), axis=None)
    #print(len(a),len(b))
    return a*d


def vectorStretch(array: np.ndarray, new_len: int):
    la = len(array)
    return np.interp(np.linspace(0, la - 1, num=new_len), np.arange(la), array)



def generate_fasta(names, sequences, output_file):
    if len(names) != len(sequences):
        raise ValueError("The number of names and sequences must be the same.")

    with open(output_file, "w") as fasta_file:
        for name, sequence in zip(names, sequences):
            fasta_file.write(f">{name}\n")
            fasta_file.write(f"{sequence}\n")




def write_fasta(sequences, names, filename):
    with open(filename, 'w') as f:
        for name, seq in zip(names, sequences):
            f.write(f'>{name}\n{seq}\n')

def allMutations(seq):
    return [seq[:i] + char + seq[i+1:] for i in range(len(seq)) for char in aminoList]


def generate_strings(original_string, substitution_chars):
    modified_sequences = [original_string[:i] + new_char + original_string[i+1:] 
                          for i, original_char in enumerate(original_string) 
                          for new_char in substitution_chars 
                          if new_char != original_char]
    
    old_characters = [original_char 
                      for i, original_char in enumerate(original_string) 
                      for new_char in substitution_chars 
                      if new_char != original_char]
    
    positions = [i 
                 for i, original_char in enumerate(original_string) 
                 for new_char in substitution_chars 
                 if new_char != original_char]
    
    new_characters = [new_char 
                      for i, original_char in enumerate(original_string) 
                      for new_char in substitution_chars 
                      if new_char != original_char]
    
    return modified_sequences, old_characters, positions, new_characters
