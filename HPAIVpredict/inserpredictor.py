
# load modules
import pandas as pd
import numpy as np
import INSERpredict3 as INSER
import argparse
from ViennaRNA import RNA

# In[18]:


#get the STDIN entry of the fasta file, 

parser = argparse.ArgumentParser(description="Process sequences in fasta for DDG analysis")

parser.add_argument("fasta", help="Name of your fasta file")
parser.add_argument("value1", type=int, default=0, help="inferior limit of window of analysis")
parser.add_argument("value2", type=int, default=50, help="superior limit of window of analysis")
parser.add_argument("value3", type=int, default=1, help="value of the start position")


args = parser.parse_args()

print("FASTA:", args.fasta)
print("inferior limit:", args.value1)
print("superior limit:", args.value2)
print("position 1:", args.value3)

samples = []
sequences = []
#transform fasta to dataframe
with open(args.fasta) as f:
    current_header = None
    current_seq = []

    for line in f:
        line = line.strip()
        
        if line.startswith(">"):
            if current_header is not None:
                samples.append(current_header)
                sequences.append("".join(current_seq))
            
            current_header = line[1:]  # quitar ">"
            current_seq = []
        else:
            current_seq.append(line)
    
    # guardar último registro
    if current_header is not None:
        samples.append(current_header)
        sequences.append("".join(current_seq))

references = pd.DataFrame({
    "sample": samples,
    "sequence": sequences
})

#print(references.head())

# get the names of the viruses

templates = references["sample"].unique()

print("your samples are: ", templates)
#mapped mutations



#launch the predictior and get the predictions  


prediction, predictioncrna = INSER.rawpredictor(references, args.value1, args.value2, args.value3)



prediction=INSER.pathway(prediction)
predictioncrna=INSER.pathway(predictioncrna)
#to save modifications in excel table
prediction.to_excel("predictiontest.xlsx")
predictioncrna.to_excel("predictioncrnatest.xlsx")

print("DDG calculation finish, excel names are predictiontest and predictioncrnatest.xlsx, please refer to R script for prediction and visualization")


###defaults for this part 

#26, 85, 986



























