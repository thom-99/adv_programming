
#parses a generic file
class Parser():
    def __init__(self,path):
      self.path = path
      self.data = []
      self.parsed = False #flag that confirmes that the file is parsed

    #parses the file, accessing it from the specified path
    def parse(self):
      try:
        with open(self.path,'r') as file:
          self.data = file.readlines()
          self.parsed = True
        print(f'{self.path} successfully parsed!')
      except FileNotFoundError:
        print(f'{self.path} was not found')

    #returns the data in the file, if the file is parsed
    def get(self):
      if self.parsed:
        return self.data
      else:
        print('you have to parse the file first using .parse method!')

    #prints the data of the file, if the file is parsed
    def show(self):
      if self.parsed:
        for line in self.data:
          print(line)
      else:
        print('you have to parse the file first using .parse method!')

    #shows the first file of the file
    def top(self):
      if self.parsed:
        print(self.data[0])


class Fasta(Parser):

    def __init__(self,path):
        super().__init__(path)
        self.processed = False

    #build a dictonary key=header value=sequence
    def process(self):
      dictonary = {}
      data = self.get()
      header = None
      sequence = ''

      for line in data:
          if line.startswith('>'):
              if header is not None:
                  dictonary[header] = sequence  # Store previous header and sequence
              header = line.strip()  # Set the new header, remove whitespaces
              sequence = ''  # Reset sequence for the new header
          else:
              sequence += line.strip()  # Add line to the current sequence (removes extra spaces/newlines)

      # Add the last header and sequence to the dictionary
      if header is not None:
          dictonary[header] = sequence
      return dictonary


    #calculates the number of distinct sequences present in the file
    def number_of_sequences(self):
        dictonary = self.process()
        return len(dictonary)

    #creates a dataframe from the fasta file
    def toDataFrame(self):
        import pandas as pd
        return pd.DataFrame(data=self.process(), index=self.process().keys())
    


  #the analyzer is an object which collects a plethora of methods to analyze a genomic sequence
#ideally the analyzer is initialized with a genomic sequence of interest
import pandas as pd


#sequence motif object can represent short DNA sequences that can be analyzed when the analyzer methods are called on a genomic sequence
#instances of this class are passed to methods to analyze the distribution of those motifs on genomic sequences

class Sequence_motif:
  def __init__(self, sequence):
    self.sequence = sequence.upper()

  #takes in input a sequence, finds distribution of motif in that sequence, useful to interact with analyzer class
  def motif_distribution(self, DNA: str):
        DNA = DNA.upper()
        list_positions = []
        index = 0
        while True:
            index = DNA.find(self.sequence, index)
            if index == -1:
                break
            list_positions.append(index)
            index += 1
        return [int(list_positions[i]+1) for i in range(len(list_positions))]


class Analyzer:

  def __init__(self, sequence):
    self.sequence = sequence.upper()

  #calculates the % of GC
  @staticmethod
  def GC_content(seq):
    GC_count = seq.count('G') + seq.count('C')
    return (GC_count/len(seq))*100

  #reverse complements the genomic sequence
  @staticmethod
  def reverse_complement(seq):
    reversecomp = ''
    for i in seq:
      if i == 'A':
        seq += 'T'
      if i == 'C':
        seq += 'G'
      if i == 'G':
        seq += 'C'
      if i == 'T':
        seq += 'A'

    return reversecomp[::-1]

  #aligns the two sequences with needleman and wunsch and returns the alignment
  #NOTE, non vi preoccupate per match_score, mismatch_score etc, teniamoli di default
  @staticmethod
  def align(seq1, seq2,gap_penalty=-2,match_score=+1,mismatch_score=-1):
    import numpy as np
    seq2 = seq2.upper()
    m = len(seq1)
    n = len(seq2)
    score_matrix = np.zeros((m+1,n+1),dtype=int)
    traceback = np.zeros((m+1,n+1),dtype=str)
    traceback[0][0] = None
    #initialization of the matrix
    for i in range(1,m+1):
      score_matrix[i][0] = i*gap_penalty
      traceback[i][0] = "up"
    for j in range(1,n+1):
      score_matrix[0][j] = j*gap_penalty
      traceback[0][j] = "down"
      #fill the matrix
    for i in range(1,m+1):
      for j in range(1,n+1):
        if seq1[i-1] == seq2[j-1]:
          score_diag = score_matrix[i-1][j-1] + match_score
        else:
          score_diag = score_matrix[i-1][j-1] + mismatch_score       #score for diagonal movement
        score_left = score_matrix[i][j-1] + gap_penalty      #changes only the column not the row
        score_up = score_matrix[i-1][j]  + gap_penalty     #changes only the row not the column
        scores = (score_diag, score_left, score_up)
        score_matrix[i][j] = max(scores)

        if score_matrix[i][j] == score_diag:
          traceback[i][j] = "diag"
        if score_matrix[i][j] == score_up:
          traceback[i][j] = "up"
        if score_matrix[i][j] == score_left:
          traceback[i][j] = "left"

    alignment_score = score_matrix[m][n]
    #now we need to compute the alignment using the traceback matrix
    aligned_seq1 = ""
    aligned_seq2 = ""
    i = len(seq1)
    j = len(seq2)
    gap = "-"
    while i>0 and j>0 :
      if traceback[i][j] == "d":
        aligned_seq1 = seq1[i-1] + aligned_seq1          #diagonal movement, added a character to both sequences
        aligned_seq2 = seq2[j-1] + aligned_seq2
        i=i-1
        j=j-1
      elif traceback[i][j] == "l":                           #movement on the left, gap on rows sequence, added a character only for columns sequence
        aligned_seq1 = gap + aligned_seq1
        aligned_seq2 = seq2[j-1] + aligned_seq2
        j = j-1
      elif traceback[i][j] == "u":
        aligned_seq1 = seq1[i-1] + aligned_seq1
        aligned_seq2 = gap + aligned_seq2
        i=i-1
    print("The score of the alignment is:", alignment_score )
    print(aligned_seq1)
    print(aligned_seq2)
    return aligned_seq1, aligned_seq2


  def motifs_analysis(self, *args : Sequence_motif):          #given a sequence and an arbitrary set of Sequence_motif instances creates a dataFrame containing the motif and some important properties
    motifs_dict = {"Motif":[],"Count":[],"Starting position": []}
    for Sequence_motif in args:
      if Sequence_motif.sequence not in motifs_dict["Motif"] and Sequence_motif.sequence in self.sequence:
        motifs_dict["Motif"].append(Sequence_motif.sequence)
        motifs_dict["Count"].append(self.sequence.count(Sequence_motif.sequence))
        motifs_dict["Starting position"].append(Sequence_motif.motif_distribution(self.sequence))   #starting positions of motifs in the DNA sequence passed to analyzer class
      else:
        raise ValueError("the motif:", Sequence_motif.sequence,"is not in the analyzed DNA sequence")
    motifs_dataframe = pd.DataFrame(motifs_dict, columns= ["Motif","Count","Starting position"])
    return motifs_dataframe

  @staticmethod
  def sequence_length(seq):
    return len(seq)

  def transcription(self):   #transcribes a DNA sequence into an RNA one
    transcribed_seq = ""
    for i in range(len(self.sequence)):
      if self.sequence[i] in ["A","C","G"]:
        transcribed_seq += self.sequence[i]
      if self.sequence[i]=="T":
        transcribed_seq += "U"
    return transcribed_seq

  def extract_subsequence(self, subsequence):   #extracts a given subsequence from the DNA sequence we are analyzing
    if subsequence not in self.sequence:
      raise ValueError("The subsequence:", subsequence, "was not found in the analyzed DNA sequence")
    else:
      print("extracting the subsequence from DNA...")
      return subsequence
    

  