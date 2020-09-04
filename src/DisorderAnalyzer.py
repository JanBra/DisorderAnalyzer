#Standardlibrary
import math
import re
import argparse
import gzip
import os
import json
import pickle
import sys
import random
from pathlib import Path

#Dependencies
import requests
import numpy as np
import pandas as pd
import scipy.stats






##############################################################      GLOBALS
REGIONMINLEN = 5

##############################################################      HELPER FUNCTIONS

def prepareMobiDBJSON():
  #Füge Komma in jede Zeile ein, um json-konform zu werden; Füge am Anfang eine Zeile mit '[' ein und am Ende eine Zeile mit ']'
  with open('disorder_UP000005640.mjson','r') as istr:
    with open('disorder_UP000005640.json','w+') as ostr:
        for i,line in enumerate(istr):
            if (i==0):
              line = "[" + line
            line = line.rstrip('\n') + ','
            print(line, file=ostr)
  #Entferne die Komma von den letzten drei Zeilen (letzter Eintrag aus Liste, Listenende("]"), und leere Zeile)
  with open('disorder_UP000005640.json','r') as myFile:
    txt = myFile.readlines()
  txt[-1] = txt[-1][0:-2] + '\n'
  txt.append("]\n")
  with open('disorder_UP000005640.json','w') as myFile:
    myFile.writelines(txt)
  os.remove('disorder_UP000005640.mjson')




#Modes: 0 ~> Bias zu Ordered; 1 ~> Bias zu Disordered; 2 ~> Probabilistisch
def liesIn(regionStart,regionEnd,peptideStart,peptideEnd,mode):
  mode = int(mode)
  if mode == 1:
    #liegt ganz in Region
    #startet vor Region & endet darin
    #startet in Region und geht darüber hinaus
    #umfasst ganze Region (startet vor und endet hinter)
    if (peptideStart >= regionStart and peptideEnd <= regionEnd) or \
        (peptideStart <= regionStart and peptideEnd <= regionEnd and peptideEnd > regionStart) or \
        (peptideStart >= regionStart and peptideEnd >= regionEnd and peptideStart <= regionEnd) or \
        (peptideStart <= regionStart and peptideEnd >= regionEnd):
        return True
    else:
        return False
  elif mode == 0:
    #liegt ganz in Region
    if (peptideStart >= regionStart and peptideEnd <= regionEnd):
      return True
    else:
      return False
  elif mode == 2:
    #Check first if it lies in Region at all(like first case for mode=1)
    if (peptideStart >= regionStart and peptideEnd <= regionEnd) or \
        (peptideStart <= regionStart and peptideEnd <= regionEnd and peptideEnd > regionStart) or \
        (peptideStart >= regionStart and peptideEnd >= regionEnd and peptideStart <= regionEnd) or \
        (peptideStart <= regionStart and peptideEnd >= regionEnd):
      countAAbeforeStart = max(regionStart - peptideStart,0)
      countAAafterEnd = max(peptideEnd - regionEnd,0)
      totalPeptideLength = peptideEnd - peptideStart
      percentageOutside = (countAAbeforeStart + countAAafterEnd)/totalPeptideLength
      outcome = random.choices([0,1],weights=[percentageOutside,1-percentageOutside])[0]
      if outcome == 0:
        return False
      else:
        return True
    else:
      return False


      

# liefert die Schnittpositionen von Trypsin aus angegebener Aminosäuresequenz
# Arguments: protein : string of protein amino acid sequence
def trypticDigestion(protein,minLen):
  #minLen = int(minLen)
  regEx = "(?<=[KR])(?!P)"
  # MITHILFE DER REGULAR EXPRESSION DIE CUT POSITIONEN BESTIMMEN
  cuts = [(match.start()) for match in re.finditer(regEx,protein)]
  cuts.insert(0,0)
  cuts.append(len(protein))
  intervals = []
  for i in range(len(cuts)-1):
    x = cuts[i]
    y = cuts[i+1]
    #intervals.append((x,y))
    if(y-x >= minLen):
      intervals.append((x,y))
  return intervals

#Merge die Intervalle, wobei die Mindestlänge 5 beträgt
# Parameter: Liste der zu mergenden Intervalle (in einem Protein)
def mergeIntervals(intervalList):
    intervalList.sort(key=(lambda x: x[0])) #Sortiere nach Start
    stack = []
    start = 0
    try:
      while intervalList[start][1] - intervalList[start][0] < 5:
        start += 1
      stack.append(intervalList[start])
    except:
      pass  #Es gab kein Start-Wert vom Intervall
    if start + 1 < len(intervalList): #Für den Fall das keine Region Länge 5 hatte
      for interval in intervalList[start+1:len(intervalList)]:
          if interval[0] <= stack[len(stack)-1][1]:
            stack[len(stack)-1][1] = interval[1]
          else:
            if interval[1] - interval[0] >= REGIONMINLEN:
              stack.append(interval)
    return stack

#Input: 
#       quality: Wertebereich: 0,1,2 ~> 0 := nur db; 1 := db+derived; 2 := db+derived+predictor
def getRegionList(protein,quality):
  regionList = []
  quality = int(quality)
  proteinFeatures = protein['mobidb_consensus']['disorder'] #spare Tipparbeit
  if "db" in proteinFeatures.keys():
    for RegionObject in proteinFeatures["db"]: #Liste von Objekten mit Methode + Liste von Regionen wird durchiteriert
      for region in RegionObject["regions"]:
        regionListItem = []
        regionListItem.append(region[0])
        regionListItem.append(region[1])
        if region[2] == 'D':
          regionList.append(regionListItem)
  if quality >= 1 and "derived" in proteinFeatures.keys():
    for RegionObject in proteinFeatures["derived"]:
      if RegionObject['method'] == 'full':
        for region in RegionObject['regions']:
          regionListItem = []
          regionListItem.append(region[0])
          regionListItem.append(region[1])
          if region[2] == 'D':
            regionList.append(regionListItem)
  if quality == 2 and "predictors" in proteinFeatures.keys():
    for RegionObject in proteinFeatures["predictors"]:
      if RegionObject["method"] == 'simple':
        for region in RegionObject['regions']:
          regionListItem = []
          regionListItem.append(region[0])
          regionListItem.append(region[1])
          if region[2] == 'D':
            regionList.append(regionListItem)
  return regionList


def getSequence(completeSequence,start,end):
  return completeSequence[start:end]

def createPeptideDataframeColumn(row):
  peptideList = []
  for peptideCut in row['TrypticPeptides']:
    peptideList.append(getSequence(row['Sequence'],peptideCut[0],peptideCut[1]))
  return peptideList

def constructDataStructure(quality,mode,minLen):
  #Öffne beide DBs;
  #Tryptisch verdauen(für jedes Protein speichere die Cut-Positionen);
  #Merge die Ungeordneten Regionen ;
  #Untersuche jedes Peptid aus den Proteinen, wenn es in einer der gemergten Regionen liegt, speicher es dementsprechend ab;
  #Guck, ob es Proteine in Uniprot gibt, die nicht in MobiDB vorkommen und nehme an, dass diese geordnet sind -> Füge davon Peptide in Datenstruktur ein
  constructNew = True
  if os.path.isfile(os.path.join(os.path.dirname(__file__),"datastructure.pkl")):
    with open(os.path.join(os.path.dirname(__file__),"datastructure.pkl"), "rb") as f:
      peptideDataStructure = pickle.load(f)
    if peptideDataStructure["Metadata"]["Quality"] == quality and peptideDataStructure["Metadata"]["Mode"] == mode and peptideDataStructure["Metadata"]["minLen"] == minLen: #Construct new datastructure if the existing file is not appropriate
      constructNew = False
  if constructNew or not os.path.isfile(os.path.join(os.path.dirname(__file__),"datastructure.pkl")):
    peptideDataStructure = {"Metadata":{}, "Peptides":{}}
    disorderCount = 0 #Zähle ungeordnete Peptide aus MobiDB-Daten
    proteinCount = 0
    orderedPeptidesCount = 0
    MobiDBAccessionSet = set() #Accessions von Proteinen aus MobiDB
    try:
      with open(os.path.join(os.path.dirname(__file__),'disorder_UP000005640.json'),'r') as mobiDBFile:
        disorderData = json.load(mobiDBFile)
      #with open('disorderTest.json','r') as mobiDBFile:
      #  disorderData = json.load(mobiDBFile)
    except:
      print("The MobiDB file could not be found. Maybe try running this program with the update option first.")
      sys.exit()
    for protein in disorderData:
      proteinCount += 1
      proteinAccession = protein['acc']
      MobiDBAccessionSet.add(proteinAccession)
      trypticCuts = trypticDigestion(protein['sequence'],minLen)
      
      #Get a flat list of the different annotated regions
      disorderedRegions = getRegionList(protein,quality)
      mergedRegions = mergeIntervals(disorderedRegions)

      disorderedPeptideIndexes = set()
      completeSetOfPeptideIndexes = set(range(len(trypticCuts)))
      for region in mergedRegions:
        for idx,peptide in enumerate(trypticCuts):  #Use Index to keep track of disordered peptides; 
          if (liesIn(region[0],region[1],peptide[0],peptide[1],mode) and idx not in disorderedPeptideIndexes):
            disorderedPeptideIndexes.add(idx)
            #disorderCount += 1
            ### ADDE PEPTID MIT DISORDER TAG
            peptideSequence = getSequence(protein['sequence'],peptide[0],peptide[1])
            if (peptideSequence not in peptideDataStructure['Peptides'].keys()):
              peptideDataStructure['Peptides'][peptideSequence] = {"Disordered": set([proteinAccession]),"Ordered": set()}
            else: #Peptid schon einmal in Datenstrukur
                peptideDataStructure['Peptides'][peptideSequence]["Disordered"].add(proteinAccession)
          if(peptide[0] > region[1]): #Peptid-Start liegt hinter Region-Ende ~> Breche ab und fahre mit nächster Region fort
            break
      orderedPeptideIndexes = completeSetOfPeptideIndexes.difference(disorderedPeptideIndexes)
      #disorderCount += len(disorderedPeptideIndexes)
      #orderedPeptidesCount += len(orderedPeptideIndexes)
      for orderedPeptide in orderedPeptideIndexes:
        # Füge geordnete Peptide mit den Indexen hinzu
        peptideCuts = trypticCuts[orderedPeptide] #Hole Schnittkoordinaten
        peptideSequence = getSequence(protein['sequence'],peptideCuts[0],peptideCuts[1])

        if (peptideSequence not in peptideDataStructure['Peptides'].keys()):
          peptideDataStructure['Peptides'][peptideSequence] = {"Disordered": set(),"Ordered": set([proteinAccession])}
          #orderedPeptidesCount += 1
        else:
          #if len(peptideDataStructure['Peptides'][peptideSequence]['Ordered'])==0: #First occurence in ordered Context of this peptide
            #orderedPeptidesCount += 1
          peptideDataStructure['Peptides'][peptideSequence]['Ordered'].add(proteinAccession)
           

    #print(peptideDataStructure)
    #raise Exception
    #AUS UNIPROT DATEN EINFÜGEN
    try:
      uniprotTable = pd.read_csv(os.path.join(os.path.dirname(__file__),"UniprotHumanProteome.tsv"),sep="\t")
    except:
      print("The Uniprot file could not be found. Maybe try running this program with the update option first.")
      sys.exit()
    #Filtern nach nicht vorgekommenen Proteinen, bevor diese Tryptisch verdaut werden
    uniprotTable = uniprotTable[~uniprotTable['Entry'].isin(MobiDBAccessionSet)]
    uniprotTable['TrypticPeptides'] = uniprotTable['Sequence'].apply(lambda x: trypticDigestion(x,minLen))
    uniprotTable['TrypticPeptides'] = uniprotTable.apply(createPeptideDataframeColumn,axis=1)

    for protein in uniprotTable.itertuples():
      
      proteinCount += 1
      TrypticPeptides = getattr(protein,"TrypticPeptides")

      #orderedPeptidesCount += len(TrypticPeptides)
      proteinAccession = getattr(protein,"Entry")
      for peptide in TrypticPeptides:
        #Füge Peptid zu Datenstruktur hinzu
        if (peptide not in peptideDataStructure['Peptides'].keys()):
          peptideDataStructure['Peptides'][peptide] = {"Disordered": set(), "Ordered":set([proteinAccession])}
          #orderedPeptidesCount += 1
        else:
          #if len(peptideDataStructure['Peptides'][peptideSequence]['Ordered'])==0:  #First occurence in ordered Context of this peptide
            #orderedPeptidesCount += 1
          peptideDataStructure['Peptides'][peptide]["Ordered"].add(proteinAccession)


    #Bestimme die Gesamtanzahl von Disordered und Ordered Peptiden (sodass jedes Peptid auch nur einmal gezählt wird)
    #for _,peptide in peptideDataStructure["Peptides"].items():
    #  #Entweder Disordered oder Ordered Einträge leer
    #  if peptide["Ordered"] and not peptide["Disordered"]:
    #    orderedPeptidesCount += 1
    #    disorderCount += 1
    #  else:
    #  elif peptide["Disordered"] and not peptide["Ordered"]:
    #    #Handle Ambiguities with Mode-Setting
    #    if mode == "0": #Bias toward order
    #      orderedPeptidesCount += 1 #Order ist vorhanden(ansonsten obere Ifs ausgelöst)
    #    elif mode == "1": #Bias toward disorder
    #      disorderCount += 1
    #    elif mode == "2":
    #      numDis = len(peptide["Disordered"])
    #      numOrd = len(peptide["Ordered"])
    #      disorderRatio = numDis / (numDis + numOrd)
    #      res = random.choices([-1,1],weights=[1-disorderRatio,disorderRatio])[0]
    #      if res == -1:
    #        orderedPeptidesCount += 1
    #      else:
    #        disorderCount += 1
    #    else:
    #      print("Mode setting unknown")

    for _,peptide in peptideDataStructure["Peptides"].items():
      orderedPeptidesCount += len(peptide['Ordered'])
      disorderCount += len(peptide['Disordered'])
        
    #Füge Metadaten zu der Datenstruktur hinzu:
    peptideDataStructure["Metadata"]["DisorderedPeptideCount"] = disorderCount
    peptideDataStructure["Metadata"]["OrderedPeptideCount"] = orderedPeptidesCount
    peptideDataStructure["Metadata"]["PeptideTotal"] = disorderCount + orderedPeptidesCount
    peptideDataStructure["Metadata"]["Quality"] = quality
    peptideDataStructure["Metadata"]["ProteinCount"] = proteinCount
    peptideDataStructure["Metadata"]["Mode"] = mode
    peptideDataStructure["Metadata"]["minLen"] = minLen

    #Save datastructure to file
    with open(os.path.join(os.path.dirname(__file__),"datastructure.pkl"), "wb+") as f:
      pickle.dump(peptideDataStructure, f, pickle.HIGHEST_PROTOCOL)

  return peptideDataStructure

# Get the peptide data from the tsv file into a usable form
# Input: String of sample peptide data
# Output: Dataframe with the Sequence and Accession Data as columns
def parsePeptideTsv(peptideFile):
  df = pd.read_csv(peptideFile[0],sep="\t")
  try:
    parsedResults = df[['Sequence','Accessions','Missed Cleavages']] #Pia Export mit Missed Cleavages Wert
    parsedResults = parsedResults.loc[parsedResults['Missed Cleavages'] == 0]
    parsedResults = parsedResults[['Sequence','Accessions']]
  except:
    parsedResults = df[['Sequence','Accessions']] #TextExporter
    parsedResults = parsedResults.drop_duplicates()
  return parsedResults

def calculateTestStatistic(ResultTable,statistic):
  if (statistic == "c2"):
    #Get ContingencyTable from ResultTable
    contigTable = [[0 for i in range(2)] for j in range(2)]
    contigTable[0][0] = ResultTable[0][0]                       #Disordered & Sample
    contigTable[0][1] = ResultTable[0][1] - ResultTable[0][0]   #Disordered & NOT_Sample
    contigTable[1][0] = ResultTable[1][0]                       #Ordered & sample
    contigTable[1][1] = ResultTable[1][1] - ResultTable[1][0]   #Ordered & NOT_Sample


    _,pValue,dof,_ = scipy.stats.chi2_contingency(contigTable,correction=False)
    assert dof == 1
  elif (statistic == "fe"):
    _,pValue = scipy.stats.fisher_exact(ResultTable)
  return pValue

#DEPRECATED DEBUG REMOVE print
def lookUp(row,datastructure,mode):
  peptideSequence = row["Sequence"]
  proteinSet = row["Accessions"] #Erst nur ein String von Liste
  proteinSet = set(proteinSet.strip('[]').replace(' ', '').split(',')) #String wird zu Liste und dann zu set konvertiert
  if peptideSequence in datastructure['Peptides'].keys():
    if proteinSet.issubset(datastructure['Peptides'][peptideSequence]['Disordered']):
      return 1
    elif proteinSet.issubset(datastructure['Peptides'][peptideSequence]['Ordered']):
      return -1
    else:
      #Handle Ambiguity with different methods
      if mode == "0":
        #PESSIMISTIC ~> Biased towards ordered peptides
        if any(protein in datastructure['Peptides'][peptideSequence]['Disordered'] for protein in proteinSet) and any(protein in datastructure['Peptides'][peptideSequence]['Ordered'] for protein in proteinSet):
          return -1
        else:
          return 0
      elif mode == "1":
        #OPTIMISTIC ~> Biased towards disordered peptides
        if any(protein in datastructure['Peptides'][peptideSequence]['Disordered'] for protein in proteinSet) and any(protein in datastructure['Peptides'][peptideSequence]['Ordered'] for protein in proteinSet):
          return 1
        else:
          return 0
      elif mode == "2":
        #ZUFALLSEXPERIMENT ~> Mache eine Anzahl von Versuchen, bei denen mit der Wahrscheinlichkeit vom Anteil von Disordered Proteinen im Set das Vorkommen als Disordered gezählt wird
        numberDisorderedProteins = 0
        numberOrderedProteins = 0
        for protein in proteinSet:
          if protein in datastructure['Peptides'][peptideSequence]['Disordered']:
            numberDisorderedProteins += 1
          if protein in datastructure['Peptides'][peptideSequence]['Ordered']:
             numberOrderedProteins += 1
        if numberDisorderedProteins == 0 and numberOrderedProteins == 0: #Keine Proteine aus Liste in Datenstruktur
          return 0
        disorderRatio = numberDisorderedProteins/(numberDisorderedProteins + numberOrderedProteins)
        return random.choices([-1,1],weights=[1-disorderRatio,disorderRatio])[0]
  else:
    return 0

def numCount(row,datastructure):
  numDis = 0
  numOrd = 0
  peptideSequence = row["Sequence"]
  proteinSet = row["Accessions"] #Erst nur ein String von Liste
  proteinSet = set(proteinSet.strip('[]').replace(' ', '').split(',')) #String wird zu Liste und dann zu set konvertiert
  if peptideSequence in datastructure['Peptides'].keys():
    for protein in proteinSet:
      if protein in datastructure['Peptides'][peptideSequence]['Disordered']:
        numDis += 1
      if protein in datastructure['Peptides'][peptideSequence]['Ordered']:
        numOrd += 1
  return (numDis,numOrd)


#Zählt wie viele geordnete und ungeordnete Peptide in der Probe vorkamen und wie die Verhältnisse in der Datenbank sind (HIER OPTIMIERBAR,INDEM DATENBANK ERGEBNISSE GESPEICHERT WERDEN???)
# Output: resultTable:              | Sample  | Database
#                       ----------------------------------
#                       Disordered  |         |
#                       ----------------------------------
#                       Ordered     |         |     
def countPeptides(sampleData, disorderDataStructure,mode):
  resultTable = [[0 for i in range(2)] for j in range(2)] #Create 2d Array
  #Nutze die gespeicherten Daten in der Datenstruktur, um die Datenbankeinträge direkt einzulesen
  resultTable[0][1] = disorderDataStructure["Metadata"]["DisorderedPeptideCount"]
  resultTable[1][1] = disorderDataStructure["Metadata"]["OrderedPeptideCount"]
  #Bestimme die Anzahlen aus der Probe
  #Vorgang: Suche in Datenstruktur nach Peptidsequenz("Sequence" aus Dataframe). Anschließend vergleiche von diesem Eintrag, ob ein Protein aus dem Accession-Eintrag aus der Probe in dem Disordered und Ordered Set vorkommen
  #Wenn alle Proteine nur in einem der Sets vorhanden sind: Erhöhe entsprechenden Counter. Ansonsten gehe vor, wie Ambiguities behandelt werden sollen.
  #results = sampleData.apply(lambda x: lookUp(x,disorderDataStructure,mode),axis=1)
  results = sampleData.apply(lambda x: numCount(x,disorderDataStructure),axis=1)

  results = pd.DataFrame(results.tolist(),index=results.index,columns=['Dis','Ord'])
  resultTable[0][0] = results['Dis'].sum()
  resultTable[1][0] = results['Ord'].sum()
  return resultTable

##############################################################      MAIN FUNCTIONS USED BY COMMAND-LINE-INTERFACE
def doUpdate():
  ###   Uniprot Data  ###
  uniprotURL = "https://www.uniprot.org/uniprot/?query=proteome:UP000005640&format=tab&columns=id,length,sequence"
  r = requests.get(uniprotURL)
  responseBody = r.text
  with open(os.path.join(os.path.dirname(__file__),"UniprotHumanProteome.tsv"),'w+') as myFile:
    myFile.write(responseBody)
  ###   MobiDB Data   ###
  mobidbURL = "https://mobidb.org/mobidb3_datasets/latest/disorder_UP000005640.mjson.gz"
  r = requests.get(mobidbURL)
  results = gzip.decompress(r.content)
  with open(os.path.join(os.path.dirname(__file__),"disorder_UP000005640.mjson"),'wb+') as myFile:
    myFile.write(results)
  prepareMobiDBJSON()
  if os.path.isfile(os.path.join(os.path.dirname(__file__),"datastructure.pkl")):
    os.remove(os.path.join(os.path.dirname(__file__),"datastructure.pkl"))
  


def doAnalyze(dataFile,statistic,quality,mode,minLen):
  ##TODO CHECK WHEN LAST RUN & WITH WHICH QUALITY(DONT'T CONSTRUCTDATASTRUCTURE WITH EVERY RUN)
  #Steps: Parse file, Count peptides(construct data structure from downloaded files with specified quality in mind), calculate the test results
  parsedPeptideContent = parsePeptideTsv(dataFile)
  peptideDataStructure = constructDataStructure(quality,mode,minLen)

  ResultTable = countPeptides(parsedPeptideContent,peptideDataStructure,mode)

  disorderedPeptidesSample = ResultTable[0][0]
  orderedPeptidesSample = ResultTable[1][0]
  disorderedPeptidesDB = ResultTable[0][1]
  orderedPeptidesDB = ResultTable[1][1]

  print("The Contingency Table looks like this:")
  print(f"\t\t|\tSample\t|\tNot found\t|\tDatabase")
  print("-----------------------------------------------------------------")
  print(f"Disordered\t|\t{disorderedPeptidesSample}\t|\t{disorderedPeptidesDB-disorderedPeptidesSample}\t\t|\t{disorderedPeptidesDB}")
  print("-----------------------------------------------------------------")
  print(f"Ordered\t\t|\t{orderedPeptidesSample}\t|\t{orderedPeptidesDB-orderedPeptidesSample}\t\t|\t{orderedPeptidesDB}")
  print("-----------------------------------------------------------------")
  print(f"Sum\t\t|\t{orderedPeptidesSample + disorderedPeptidesSample}\t|\t{orderedPeptidesDB-orderedPeptidesSample+disorderedPeptidesDB-disorderedPeptidesSample}\t\t|\t{orderedPeptidesDB+disorderedPeptidesDB}")


  pValue = calculateTestStatistic(ResultTable,statistic)
  print(f"The p-Value is: {pValue}")
  print(f"If your confidence level is higher than this value, the null hypothesis can be rejected and disordered peptides were over- or underrepresented.")
  confidenceLevel = float(input("Please enter your confidence level: "))
  if(pValue < confidenceLevel):
    print("Null hypothesis can be rejected. Disordered peptides were significantly", end=" ")
    if (disorderedPeptidesSample/(disorderedPeptidesSample + orderedPeptidesSample) > disorderedPeptidesDB/(disorderedPeptidesDB + orderedPeptidesDB)):
      print("overrepresented.")
    else:
      print("underrepresented.")
  else:
    print("Null hypothesis can not be rejected.")
  




##############################################################      COMMAND-LINE-INTERFACE
#########################################         MAYBE ADD CONFIDENCE LEVEL AS ARGUMENT // ARGUMENT METHOD TO RESOLVE AMBIGUITY // OPTION DONT SAVE DATASTRUCTURE
parser = argparse.ArgumentParser(description="...")
subparsers = parser.add_subparsers(title="Subcommands",dest='cmd',description="Either analyze a file with this program or update the used databases for analysis. You can use the help command again together with the subcommand to get detailed information about the Arguments needed by the subcommands.(I.e. DisorderAnalyzer.py analyze -h)")
#Optional Arguments
# Two Sub-commands: Analyze & Update

parser_analyze = subparsers.add_parser('analyze',help="Use this to analyze your tsv file")
parser_analyze.add_argument("-f","--tsvFile",nargs=1,type=Path,help="The tsv file that is supposed to be analyzed. The program was originally developed to work with the PIA peptide export. However the only requirements the file must fulfill is the containment of the following entries: 'Sequence' -> amino acid sequence of peptides in sample,'Accessions' -> results of peptide identification, the accessions of possible proteins in which the peptide occurs.")
parser_analyze.add_argument("-s","--statistic",default="c2",choices=["c2","fe"],help="Choose the statistical test that you want to use for your analysis. You can choose either Chi-Square test with 'c2' or Fisher's exact test with 'fe'. Default is Chi-Square test.")
parser_analyze.add_argument("-q","--quality",default=2,choices=["0","1","2"],help="Specify the data quality to be used for the analysis. Possible Values: 0,1 or 2. 0 ~> Only manually curated data, 1 ~> Additional to manually curated data, use data derived from experiments, 2 ~> Additionally use data generated from disorder predictors. It is strongly recommended to use 2, because the amount of data without predictors is very limited at the moment.")
parser_analyze.add_argument("-m","--mode",required=True,choices=["0","1","2"],help="This setting specifies the definition of the 'lies in' relationship between peptides and disordered regions. 0 ~> The peptide must be fully enclosed in a disordered region to be classified as disordered; 1 ~> Atleast one amino acid of a peptide must reach into a disordered region; 2 ~> probabilistic approach: the peptide is classified as disordered with the probability of the ratio of amino acids in the region vs outside the region.")
parser_analyze.add_argument("-l","--minLen",required=True,help="The minimum length used for the peptide identification.")
parser_update = subparsers.add_parser('update',help="Use this to update the databases. If you just downloaded the program you must run this command once.(No further arguments with this command)")

args = parser.parse_args()

#"Select the method to resolve ambiguity whether the peptide is sampled from an disordered or ordered protein if both may be possible. Possible values: 0: assume the protein was ordered , 1: assume the protein was disordered, 2: select randomly with the probability of counting as disordered peptide being the ratio of possible disordered proteins to disordered and ordered proteins."

if args.cmd == 'update':
  doUpdate()
elif args.cmd == 'analyze':
  minLen = max(1,int(args.minLen))
  doAnalyze(args.tsvFile,args.statistic,args.quality,args.mode,minLen)
elif args.cmd == None:
  print("Use the help command to display the commands of this program. (I.e. DisorderAnalyzer.py -h)")