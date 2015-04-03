import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

gtfLoc=sys.argv[1]
chrLoc=sys.argv[2]



class mdict(dict):

    def __setitem__(self, key, value):

        self.setdefault(key, []).append(value)


#Open database stuff

refseqToCommon=dict()
commonToRefseq=mdict()

startCodon=dict()
stopCodon=dict()
exons=mdict()
chromosome=dict()
strand=dict()
def getRevComp(seq):
	seq=seq[::-1]
	revComp=""
	for s in seq:
		if s=='G':
			revComp+='C'
		if s=='C':
			revComp+='G'
		if s=='T':
			revComp+='A'
		if s=='A':
			revComp+='T'
		if s not in ['A','C','T','G']:
			revComp.append(s)

	return revComp
for line in open(gtfLoc):
	line=line.split('\t')
	if len(line[0])<3:
		gene = line[8].split('pt_id \"')[1].split('\"; tss_')[0]
	
		if line[2]=='start_codon':
			startCodon[gene] = int(line[3])
			common= line[8].split('gene_id \"')[1].split('\"; transc')[0]
			refseqToCommon[gene]=common
			commonToRefseq[common]=gene
			chromosome[gene] = line[0]
			strand[gene]=line[6]
		if line[2]=='stop_codon':
			stopCodon[gene] = int(line[3])	
		if line[2]=='exon':
			exons[gene]=[int(line[3]),int(line[4])]
			

		
seqs=dict()
seq=''
chrNum=""
seqLen=0
cdsStart=dict() 
cdsStop=dict() 
cdsLen=dict()
mRNALen=dict()

genesToUse=[]
transcripts=mdict()

thisCDSStart=-1;thisCDSStop=-1;seqLen=0
for line in open ('transcripts.gtf'):	
	try:
		line=line.split('\t')

		if len(line[0])<6:
			if line[2]=='transcript':
				info=line[8].split(';')
				gene=info[0].split('\"')[1].replace('\'','')
				transcript=info[1].split('\"')[1].replace('\'','')
				fpkm=info[2].split('\"')[1].replace('\'','')
				if transcript[1] =='M':
					transcripts[gene]=fpkm + ',' + transcript
		

				if seqLen>0:
				
					
					
					if strand[geneToWrite]=='-': #Flip variables if strand is negative
						seq=getRevComp(seq)
						thisCDSStart=seqLen-thisCDSStart
						thisCDSStop=seqLen-thisCDSStop
					
					mRNALen[geneToWrite.split('.')[0]]=seqLen
					seqs[geneToWrite.split('.')[0]]=seq
					if strand[geneToWrite]=='+':
						thisCDSStart+=2
					cdsStart[geneToWrite.split('.')[0]]=thisCDSStart-1
					cdsStop[geneToWrite.split('.')[0]]=thisCDSStop
					mRNALen[geneToWrite.split('.')[0]]=seqLen
					cdsLen[geneToWrite.split('.')[0]]=thisCDSStop-thisCDSStart
					
					
					
					seq=""
					thisCDSStart=-1;thisCDSStop=-1;seqLen=0
				if line[0] != chrNum:
					chrSeq = SeqIO.read(open(chrLoc + line[0] + ".fa"), "fasta")
					chrNum=line[0]
					print 'Chromosome ' + chrNum + '...'
				geneToWrite=transcript
			if line[2]=='exon':
				
				if startCodon[transcript] >= int(line[3]) and startCodon[transcript] <= int(line[4]):
					thisCDSStart=seqLen + startCodon[transcript] - int(line[3])
				if stopCodon[transcript] >= int(line[3]) and stopCodon[transcript] <= int(line[4]):
					thisCDSStop=seqLen + stopCodon[transcript] - int(line[3])
				seqLen+=int(line[4]) - int(line[3])
			
				seq=seq+str(chrSeq[int(line[3]):int(line[4])].seq) 
	except KeyError: pass

writeIndex=open('database/chosenTranscripts','w')
writeDatabase=open('database/geneInfo','w')

print 'Printing database'
t=0
f=0

for gene in transcripts.keys():
	trans=transcripts[gene]
	trans.sort()
	transcript=trans[-1].split(',')[1].split('.')[0]
	
	
	t+=1
	try:
		
		thisStart=str(cdsStart[transcript])
		thisStop=str(cdsStop[transcript])
		thisSeq=seqs[transcript]
		thisLen=str(mRNALen[transcript])
		writeIndex.writelines('\n>' + transcript + '\tcommon=' + gene + '\tcdsStart=' + thisStart + '\tcdsStop=' + thisStop + '\n' +seqs[transcript])
		writeDatabase.writelines('\n' + transcript + '\t' + gene + '\t' + thisStart + '\t' + thisStop + '\t' + thisLen + '\t' )
		
	except KeyError: 
		f+=1
writeIndex.close()
writeDatabase.close()

