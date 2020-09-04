import sys
import os.path
import re

def remBracket(seq):
  newSeq = seq
  if(newSeq[0] == '.'):
    newSeq = newSeq[1:]
  
  while newSeq.find('(') != -1:
    startIdx = newSeq.find('(')
    counter = 0
    for i,char in enumerate(newSeq[startIdx:]):
      if char=='(':
        counter += 1
      elif char==')':
        counter -= 1
        if counter == 0:
          endIdx = i
          modSeq = newSeq[:startIdx] + newSeq[startIdx + endIdx + 1:]
          newSeq = modSeq

          break
  return '"' + newSeq + '"'



if __name__ == "__main__":
  os.path.join(os.path.dirname(__file__),"datastructure.pkl")
  if (len(sys.argv) > 2 or not os.path.isfile(os.path.join(os.path.dirname(__file__),sys.argv[0]))):
    print(len(sys.argv))
    raise Exception
  else:
    count = 0
    with open(sys.argv[1],'r') as istr:
      index = sys.argv[1].rfind('.')
      with open(sys.argv[1][0:index] + "_mod" + sys.argv[1][index:],'w') as ostr:
        seqIndex = 0
        accIndex = 0

        lines = istr.readlines()
        #Edit Column Header Names
        fields = lines[0].split('\t')
        
        accIndex = fields.index('accessions')
        seqIndex = fields.index('sequence')

        fields = map(lambda x: x.capitalize() , fields)

        lines[0] = f"\t".join(fields)

        newLines = []
        newLines.append(lines[0])
        #Edit the formatting of the protein accessions and sequence
        #TODO
        for line in lines[1:]:
          fields = line.split('\t')
          #Edit Sequence
          #fields[seqIndex] = '"' + re.sub(r'\([^)]*\)','',fields[seqIndex]) + '"'
          fields[seqIndex] = remBracket(fields[seqIndex])
          if(len(re.findall("(?<=[KR])(?!P)",fields[seqIndex])) > 1): #Match with missed cleavage gets ommitted
            continue
          #Edit Accessions TODO
          prots = fields[accIndex].split(';') #Split the accessions of different proteins
          prots = list(map(lambda x: x.split('|')[1],prots))
          fields[accIndex] = '"[' + ','.join(prots) + ']"'
          myLine = f"\t".join(fields)
          newLines.append(myLine)
        ostr.writelines(newLines)
      