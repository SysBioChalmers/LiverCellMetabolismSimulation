#!/usr/local/bin/python3.3
import os
os.chdir('D:/stories/Glutamine/aminoacidcomposition/proteomicsToAAComp/ens2uni')
FgeneKCat = open('result.txt', 'w')

genesForMaping = open('input.txt').read().splitlines()


bMap = {}
with open("ens2uni.txt") as f:
   for line in f:
      currentVals = line.split("\t")
      if len(currentVals)>0:
         key = currentVals[0]
         if len(currentVals)>1:
            test = currentVals[1]
            val = test.split("\n")[0]
         else:
            val = "-"
         bMap[key] = val


for line in genesForMaping:
   currentVals = line.split("\t")
   if currentVals[0] in bMap:
      if bMap[currentVals[0]]:
         FgeneKCat.write(bMap[currentVals[0]] + "\t" + str(currentVals[1]) + "\n")
FgeneKCat.close()
