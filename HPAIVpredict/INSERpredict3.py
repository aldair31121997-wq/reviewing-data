from ViennaRNA import RNA
import pandas as pd
import numpy as np
import statsmodels.api as sm
from Bio.Seq import Seq

# Function for vrna/crna analysis from the dataframe containing the sequences/sample name
# 'a' sets the start position for the analysis and 'b' sets the stop position for the analysis (we used 1030 and 1100). 'c' sets the position corrector and depends on your sample (we used 983)

def rawpredictor(references, a, b, c):

 # now getting all the variables that describe the polimerase sliding trough the RNA
 templates=references["sample"].unique()
 sisterstrand="CCCCCCCCC"
 m=1
 infi=0
 infi2=0
 decalagesplus=pd.DataFrame(columns=['interna',"DeltaG"])
 prediction= pd.DataFrame(columns=['cRNA', 'vRNA', 'imperectpair', 'decalage', 'deltaG', "pos", "insertion","sample","DDG","structure"]) 
 predictioncrna= pd.DataFrame(columns=['cRNA', 'vRNA', 'imperectpair', 'decalage', 'deltaG', "pos", "insertion","sample","DDG","structure"]) 
 size=9
 size1=10
 size2=10
 # first iterative step reads the list of templates sequence by sequence
 for v in (range(len(templates))):
  line2="" 
  rnaseq = references["sequence"][v]
  rnaseq = rnaseq[a:b]
  rnaseq = rnaseq.replace('T', 'U')
  # if you want the vRNA
  rnaseq2 = str(Seq(rnaseq).reverse_complement())
  rnaseq2 = rnaseq2.replace('T', 'U')
  # to get the complementry strand of vRNA
  rnaseqcomp = str(Seq(rnaseq).complement())
  rnaseqcomp2 = str(Seq(rnaseq2).complement())
  name = references["sample"][v]
  swindow= 1
  # This variable is to define the two interacting sequences, the & tells the algorithm that we are comparing two different sequences

############################################################################################################

  #read the length of the RNA seq
  Length = len(rnaseq)

  #This is the start point of iteration
  i = 1

  # End the replication bubble sequence and add 1, because sequence count starts at 0
  End = int((Length - size) / swindow)
 
  # Start iteration between sequences to form structure
  for i in range(9, End-9):

   decalage1=sisterstrand
   structure=["j","j","j","j","j","j","j","j","j","&","j","j","j","j","j","j","j","j","j"]
   structuredec1=["j","j","j","j","j","j","j","j","j","&","j","j","j","j","j","j","j","j","j"]
   structuredec1m=["j","j","j","j","j","j","j","j","j","&","j","j","j","j","j","j","j","j","j"]

  #here you define the catalitic site inside the RNA pol 
   cataliticsite = rnaseq[i:i+size]
   cataliticsite2 = rnaseq2[i:i+size]
   decalageminus1= rnaseq[i-1:i+size-1]
   if i>1:
    decalageminus2= rnaseq[i-2:i+size-2]
   else:
    decalageminus2=decalageminus1
   invsisterstrand = rnaseqcomp[i:i+size]
   invsisterstrand2 = rnaseqcomp2[i:i+size]
   #inverse the reversed sister strand to have the proper sense 5'3'
   sisterstrand = invsisterstrand[::-1]
   sisterstrand2 = invsisterstrand2[::-1]
   #this for creates the interaction between template and product in ((..)) notation, this is during vRNA syntesis
   for j in range(len(cataliticsite)):
     if (cataliticsite[j]=="A") & (sisterstrand[size-(j+1)] == "U"):
         structure[j]= "("
         structure[j+size+1]=")"
     elif (cataliticsite[j]=="U") & (sisterstrand[size-(j+1)] == "A"):
         structure[j]= "("
         structure[j+size+1]=")"
     elif (cataliticsite[j]=="G") & (sisterstrand[size-(j+1)] == "C"):
         structure[j]= "("
         structure[j+size+1]=")"
     elif (cataliticsite[j]=="C") & (sisterstrand[size-(j+1)] == "G"):
         structure[j]= "("
         structure[j+size+1]=")"
     else:
         structure[j]="."
         structure[j+size+1]="."
      
 ############################################################################################
   line="---------"
   perfectpair=0
   for s in range(size1):
    structurep= ["j","j","j","j","j","j","j","j","j","&","j","j","j","j","j","j","j","j","j"]
    cataliticsitep = rnaseq[i+s:i+size+s] 
    for h in range(len(cataliticsitep)):
     if (cataliticsitep[h]=="A") & (sisterstrand[size-(h+1)] == "U"):
         structurep[h]= "("
         structurep[18-h]=")"
     elif (cataliticsitep[h]=="U") & (sisterstrand[size-(h+1)] == "A"):
         structurep[h]= "("
         structurep[18-h]=")"
     elif (cataliticsitep[h]=="G") & (sisterstrand[size-(h+1)] == "C"):
         structurep[h]= "("
         structurep[18-h]=")"
     elif (cataliticsitep[h]=="C") & (sisterstrand[size-(h+1)] == "G"):
         structurep[h]= "("
         structurep[18-h]=")"
     elif (cataliticsitep[h]=="U") & (sisterstrand[size-(h+1)] == "G"):
         structurep[h]= "("
         structurep[18-h]=")"
     elif (cataliticsitep[h]=="G") & (sisterstrand[size-(h+1)] == "U"):
         structurep[h]= "("
         structurep[18-h]=")"
     else:
         structurep[h]="."
         structurep[18-h]="."
    structurep="".join(structurep)
    deltaGdecalp=RNA.eval_structure_simple(cataliticsitep+"&"+sisterstrand, structurep)
    decalagesplus.loc[s]= [structurep] + [deltaGdecalp]
    # Extract all data of predictions if the DeltaG is favorable
    prediction.loc[infi]= [cataliticsite] + [sisterstrand[::-1]] + [cataliticsitep+"&"+sisterstrand] + [s] + [decalagesplus["DeltaG"][s]] + [i+c] + [cataliticsite[0:s]] + [name] + [(decalagesplus["DeltaG"][s] - perfectpair)] + [structurep]
    infi=infi+1
    # Get the deltaG of the perfect pairing for calculation of DDG
    if s==0:
     perfectpair= deltaGdecalp
     # Analysis of the cRNA syntesis: The first external FOR iterates over each possible gap
     # The second INTERNAL FOR iterates over each nucleotide inside the catalitic site to determine the structure
   perfectpair=0
   for l in range(size2):
    structurem= ["j","j","j","j","j","j","j","j","j","&","j","j","j","j","j","j","j","j","j"]
    cataliticsitem = rnaseq2[i+l:i+size+l] 
    for h in range(len(cataliticsitem)):
     if (cataliticsitem[h]=="A") & (sisterstrand2[size-(h+1)] == "U"):
         structurem[h]= "("
         structurem[18-h]=")"
     elif (cataliticsitem[h]=="U") & (sisterstrand2[size-(h+1)] == "A"):
         structurem[h]= "("
         structurem[18-h]=")"
     elif (cataliticsitem[h]=="G") & (sisterstrand2[size-(h+1)] == "C"):
         structurem[h]= "("
         structurem[18-h]=")"
     elif (cataliticsitem[h]=="C") & (sisterstrand2[size-(h+1)] == "G"):
         structurem[h]= "("
         structurem[18-h]=")"
     elif (cataliticsitem[h]=="U") & (sisterstrand2[size-(h+1)] == "G"):
         structurem[h]= "("
         structurem[18-h]=")"
     elif (cataliticsitem[h]=="G") & (sisterstrand2[size-(h+1)] == "U"):
         structurem[h]= "("
         structurem[18-h]=")"
     else:
         structurem[h]="."
         structurem[18-h]="."
    structurem="".join(structurem)
    deltaGdecalm=RNA.eval_structure_simple(cataliticsitem+"&"+sisterstrand2, structurem)
    cataliticsiteinv= cataliticsite2[::-1]
    vrnainser=cataliticsiteinv[0:l]
    line=line+"-"
    predictioncrna.loc[infi2]= [cataliticsite] + [sisterstrand[::-1]] + [cataliticsitem+"&"+sisterstrand2] + [l+1] + [deltaGdecalm] + [i+(c-9)] + [str(Seq(cataliticsite2[0:l]).complement())] + [name] + [(deltaGdecalm - perfectpair)] + [structurem]
    infi2=infi2+1
    if l==0:
     perfectpair= deltaGdecalm
#############################################################################33
   nucleotide= rnaseq[i+8]
   line2=line2+"-"
   m = m + 1
 # correct the gap of CRNA
 predictioncrna["decalage"]=(predictioncrna["decalage"]-1)
    
 return(prediction, predictioncrna)

### Function for aditional treatment of thermodinamical pathways

def pathway(prediction):
    
 prediction["neighbor"]= 100.0

 for i in range(len(prediction)-7):
  if ((prediction["decalage"][i] > 1) & (prediction["decalage"][i] < 9)):
       prediction.loc[i,"neighbor"]= prediction["DDG"][i-1] + prediction["DDG"][i+1]
  if (prediction["decalage"][i] ==1):
       prediction.loc[i,"neighbor"]= prediction["DDG"][i+1]
  if (prediction["decalage"][i] ==9):
       prediction.loc[i,"neighbor"]= prediction["DDG"][i-1]

 # create a new column for the neighbor variable direct block
 prediction["neighborblock"]= 100.0
 prediction["equilibrated"]=1.0
 prediction["pathway"]=1.0
 prediction["pathwayb"]=1.0
 prediction["falsepos"]=prediction["pos"]
 for i in range(len(prediction)-1):      
    
  if ((prediction["decalage"][i] > 1) & (prediction["decalage"][i] < 9)):
       prediction.loc[i,"neighborblock"]= prediction["DDG"][i-1]
  if (prediction["decalage"][i] ==1):
       prediction.loc[i,"neighborblock"]= prediction["DDG"][i]
  if (prediction["decalage"][i] ==9):
       prediction.loc[i,"neighborblock"]= prediction["DDG"][i-1]
  # create a new variable to iterate over the neighborhood and search for previous stable states
  acumulation=1
  if (prediction["decalage"][i]>3):
  
   for j in range(prediction["decalage"][i]-1):
        acumulation=acumulation*prediction["deltaG"][i-(j+1)]
        # this line compares the new point with the old one and records the biggest

        prediction.loc[i,"equilibrated"]=acumulation
       # finally create the unmapping variable analysis

###########################################################"""
  neighborc=1
  neighbord=1
  if (prediction["decalage"][i]>3):
  
   for j in range(prediction["decalage"][i]):
        # this line compares the new point with the old one and records the biggest
        if (neighborc < prediction["deltaG"][i-(j+1)]):
         neighborc=(prediction["deltaG"][i-(j+1)]+0.1)
        # do the same but with DDG
        if (neighbord < prediction["DDG"][i-(j+1)]):
         neighbord=(prediction["DDG"][i-(j+1)]+0.1)

   prediction.loc[i,"pathway"] = neighborc+0.1
   prediction.loc[i,"pathwayb"] = neighbord+0.1

 # create a new column with the sequence data
 prediction['template'] = prediction['insertion'].astype(str).str[0]

 # end for creating the data regarding the template and the pseudo mapping simmulation
 for i in range(len(prediction)):
  if (i<1):
      dummypos=prediction["pos"][i]

  if ((i > 1) & (prediction["template"][i]=="A") ):
    
    if ((prediction["template"][i-2]!="A") & (prediction["decalage"][i-2]!=0)):
         dummypos=prediction["pos"][i]
    if ((prediction["template"][i-2]=="A") | (prediction["decalage"][i-2]==0)):
         prediction.loc[i,"falsepos"]=dummypos
        
  # lastly, add the match mismatch rules, with different categorical variables A, B, C, D,

 prediction["matchmismatch"]= "NOVALUE"
 
 for i in range(len(prediction)-7):
      # first case: penality mismatches D
      if ((prediction["imperectpair"][i][0]=="A") & (prediction["imperectpair"][i][18]=="G")):
       prediction.loc[i,"matchmismatch"]="D"
      if ((prediction["imperectpair"][i][0]=="G") & (prediction["imperectpair"][i][18]=="A")):
       prediction.loc[i,"matchmismatch"]="D"
      if ((prediction["imperectpair"][i][0]=="G") & (prediction["imperectpair"][i][18]=="G")):
       prediction.loc[i,"matchmismatch"]="D"
      if ((prediction["imperectpair"][i][0]=="C") & (prediction["imperectpair"][i][18]=="C")):
       prediction.loc[i,"matchmismatch"]="D"
      if ((prediction["imperectpair"][i][0]=="A") & (prediction["imperectpair"][i][18]=="A")):
       prediction.loc[i,"matchmismatch"]="D"    
      # second case: tolerated mismatches C-A, U-U, C-U VALUE B
      if ((prediction["imperectpair"][i][0]=="C") & (prediction["imperectpair"][i][18]=="A")):
       prediction.loc[i,"matchmismatch"]="C"
      if ((prediction["imperectpair"][i][0]=="A") & (prediction["imperectpair"][i][18]=="C")):
       prediction.loc[i,"matchmismatch"]="C"
      if ((prediction["imperectpair"][i][0]=="U") & (prediction["imperectpair"][i][18]=="U")):
       prediction.loc[i,"matchmismatch"]="C"
      if ((prediction["imperectpair"][i][0]=="C") & (prediction["imperectpair"][i][18]=="U")):
       prediction.loc[i,"matchmismatch"]="C"
      if ((prediction["imperectpair"][i][0]=="U") & (prediction["imperectpair"][i][18]=="C")):
       prediction.loc[i,"matchmismatch"]="C"
      # third case: matches with no plus value A-U, G-C,
      if ((prediction["imperectpair"][i][0]=="C") & (prediction["imperectpair"][i][18]=="G")):
       prediction.loc[i,"matchmismatch"]="B"
      if ((prediction["imperectpair"][i][0]=="G") & (prediction["imperectpair"][i][18]=="C")):
       prediction.loc[i,"matchmismatch"]="B"
      if ((prediction["imperectpair"][i][0]=="A") & (prediction["imperectpair"][i][18]=="U")):
       prediction.loc[i,"matchmismatch"]="B"
      if ((prediction["imperectpair"][i][0]=="U") & (prediction["imperectpair"][i][18]=="A")):
       prediction.loc[i,"matchmismatch"]="B"
      if ((prediction["imperectpair"][i][0]=="U") & (prediction["imperectpair"][i][18]=="G")):
       prediction.loc[i,"matchmismatch"]="B"
      if ((prediction["imperectpair"][i][0]=="G") & (prediction["imperectpair"][i][18]=="U")):
       prediction.loc[i,"matchmismatch"]="B"
      # fourth case: create adenine insertions inside a poli-A stretch
      if ((prediction["insertion"][i]=="A") & (prediction["insertion"][i+3]=="AAAA") ):
       prediction.loc[i,"matchmismatch"]="A"

 return(prediction)   


############# analysis for window of 10 nt ########################################################################################
