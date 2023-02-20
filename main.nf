#!/usr/bin/env nextflow 

params.cutoff = 100

process filterSeq { 
  input: 
    path input_file 
    val cutoff

  output: 
    path 'output.fasta' 

  script: 
    """
      #!/usr/bin/env python3
      from Bio import SeqIO
      data=SeqIO.parse('$input_file', 'fasta')
      sequenceList=[]
      for i in data:
        if (len(i.seq)) > $cutoff:
          sequenceList.append(i)
      SeqIO.write(sequenceList, "output.fasta", "fasta")
    """
} 

workflow { 
  inputFile = Channel.fromPath(params.input) 
  filteredData = filterSeq(inputFile, params.cutoff) 
}
