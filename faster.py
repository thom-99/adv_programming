
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
  @staticmethod
  def align(seq1, seq2, gap_penalty=-2, match_score=2, mismatch_score=-2):
    import numpy as np
  #matrix dimension

    m = len(seq1)
    n = len(seq2)
    score_matrix = np.zeros((m+1, n+1), dtype=int)
    traceback_matrix = np.zeros((m+1, n+1), dtype=str)
    max_score=0
    max_position=(0,0)

    #initialize the score matrix
    for i in range(m+1):
      score_matrix[i,0]= gap_penalty*i
      traceback_matrix[i,0] = "u"  #the matrix gives back "u" when there is a gap in seq2
    for j in range(n+1):
      score_matrix[0,j]= gap_penalty*j
      traceback_matrix[0,j]= "l" #the matrix gives back "l" when there is a gap in seq1
    traceback_matrix[0,0] = None  #start position
                                  #"u" stands for up, "l" stands for left, "d" stands for diagonal

    #fill the score matrix

    for i in range(1, m+1):
      for j in range(1, n+1):
        if seq1[i-1] == seq2[j-1]:
          diag = score_matrix[i-1][j-1] + match_score
        else:
          diag = score_matrix[i-1][j-1] + mismatch_score

        up = score_matrix[i-1][j] + gap_penalty  # gap in seq2
        left = score_matrix[i][j-1] + gap_penalty  # gap in seq1

        score_matrix[i][j]=max(0, diag, up, left)

        if score_matrix[i][j]>=max_score:
          max_score = score_matrix[i][j]
          max_position = (i,j)


        #direction of the traceback matrix
        if score_matrix[i][j]== 0:
          traceback_matrix[i][j]= None
        elif score_matrix[i][j] == diag:
          traceback_matrix[i][j] = 'd'  #match or a mismatch in the sequence
        elif score_matrix[i][j] == up:
          traceback_matrix[i][j] = 'u'    #gap in seq2
        else:
          traceback_matrix[i][j] = 'l'  #gap in seq1

    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = max_position
    while i > 0 and j > 0 and score_matrix[i][j]>0:
      if traceback_matrix[i][j] == "d":
        aligned_seq1 = seq1[i-1] + aligned_seq1
        aligned_seq2 = seq2[j-1] + aligned_seq2
        i = i-1
        j = j-1
      elif traceback_matrix[i][j] == "u":
        aligned_seq1 = seq1[i-1] + aligned_seq1
        aligned_seq2 = "-" + aligned_seq2
        i = i-1
      elif traceback_matrix[i][j] == "l":
        aligned_seq1 = "-" + aligned_seq1
        aligned_seq2 = seq2[j-1] + aligned_seq2
        j = j-1

    return aligned_seq1, aligned_seq2, score_matrix[max_position]


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
    

  