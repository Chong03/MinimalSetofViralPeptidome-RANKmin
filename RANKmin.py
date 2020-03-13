from Bio import SeqIO
import pandas as pd
from tqdm import tqdm
import os
from timeit import default_timer as timer

fileNumber = 0

# full sequence 
fullSequence = list(SeqIO.parse("dengueprotSeqPilot","fasta"))

# substring
subStrings = pd.read_csv("nr_kmers_9800", header = None)    

  
s = timer()

# repeat until substring file is empty
while(len(subStrings[0]) != 0): 
    
    oneSequence = fullSequence[fileNumber].seq
    
    occurence = 0
    
    for index, subString in tqdm(enumerate(subStrings[0])):            
        if(subString in oneSequence):
            occurence = occurence + 1
            subStrings.drop(index, inplace=True)
            
    # save fileZ        
    if(occurence != 0):
        fileZ = open("fileZ.txt","a")
        fileZ.write(str(fullSequence[fileNumber].id) + ';' + str(occurence) + '\n')
        fileZ.close()
                
    # save match to file
    fileO = open("fileO.txt","a")
    fileO.write(str(fullSequence[fileNumber].id) + ';' + str(occurence) + '\n')
    fileO.close()
    
    # increment        
    fileNumber = fileNumber + 1
    
    # save remaining substring
    subStrings.to_csv("try/outputfile_pilotDengue"+str(fileNumber), index=False, header=False)
    
    if (os.stat("try/outputfile_pilotDengue"+str(fileNumber)).st_size == 0):
        e = timer()
        time = (e-s)/60
        print("Total time used : " + str(time) + " min")
        #break
        #print("Stop reading file as it is empty")
    else:
        subStrings = pd.read_csv("try/outputfile_pilotDengue"+str(fileNumber), header=None)

###
          
highestMatching = pd.read_csv("fileZ.txt", delimiter=';', header=None, names=['id','match']).sort_values(by=['match'],ascending=False)
# read substring again
subStrings = pd.read_csv("nr_kmers_9800", header = None)    

n=len(fullSequence)
r = 0
s = timer()
# repeat until substring file is empty
while(len(subStrings[0]) != 0): 
    
    # take highest matching index
    for i in range(n):
        if fullSequence[i].id == str(highestMatching.id.iloc[r]):
            highest_index = i
            
    # take highest matching sequence
    oneSequence = fullSequence[highest_index].seq
    
    occurence = 0
        
    for index, subString in tqdm(enumerate(subStrings[0])):            
        if(subString in oneSequence):
            occurence = occurence + 1
            subStrings.drop(index, inplace=True)
    
    # save fileZ        
    if(occurence != 0):
        fileZ = open("Z.txt","a")
        fileZ.write(str(fullSequence[highest_index].id) + ';' + str(occurence) + '\n')
        fileZ.close()
    
    # save match to file
    fileO = open("O.txt","a")
    fileO.write(str(fullSequence[highest_index].id) + ';' + str(occurence) + '\n')
    fileO.close()
    
    subStrings.reset_index(drop = True, inplace=True)
    
    r = r + 1

e = timer()
time = (e-s)/60
print("Total time used : " + str(time) + " min")