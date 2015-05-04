readCutoff=3
individReadCutoff=2500
posxmin=-20
posxmax=125
ribo5Add=14

import matplotlib.pyplot as plt
import numpy as np
import random
import os

class mdict(dict):

    def __setitem__(self, key, value):
        """add the given value to the list of values for this key"""
        self.setdefault(key, []).append(value)

def makeHistoList(list,xstart=0,xstop=None,normalizeToOne=True):
	if xstop==None:
		xtop=max(list)
	list.sort();listLen=len(list)
	out=[]
	pos=0
	for i in range(xstart,xstop):
		count=list.count(i)
		
		out.append(count)
	if normalizeToOne:
		s=np.sum(out)
		newY=[]
		for i in out:
			newY.append(float(i)/s*100)
		out=newY
	return out
	
def writePosGraphs(listOfPos,xmin,xmax,name,gene="",normalize=1):

	y=makeHistoList(listOfPos,xstart=xmin,xstop=xmax,normalizeToOne=True)
	
	
	x=range(xmin,xmax)
	
	
	if len(gene)>0:
		fill=[1]*(cdsStop[gene]-cdsStart[gene])
		plt.fill_between(range(cdsStart[gene],cdsStop[gene]),0,max(y),color=(0,0,0))
		
	write=open(name + '.txt','w')
	write.writelines('Position on mRNA\tPercentage of reads')
	
	for i in range(len(y)):
		write.writelines('\n' + str(x[i]) + "\t" + str(y[i]))
		
	write.close()
	plt.plot(x,y,color=(1,0,0))
	
	plt.xlabel('Position in mRNA')
	plt.ylabel('Percentage of reads')
	plt.axis(xmin=xmin,xmax=xmax,ymax=max(y))
	plt.savefig(name + '.svg')
	plt.savefig(name + '.png')
	plt.clf()
	
	return y

def loadReads(folder,fileName):
	file = folder + '/' +fileName 
	print '\t   ' + fileName.strip('.map')
	thisCount=dict()
	c=0
	geneLevel=dict()
	individGenes=mdict()
	try:
		for line in open(file):
			line=line.split("\t")
			
			if len(line)>1:
				pos=int(line[3]) + ribo5Add
				gene=line[2]
				try:
					
					if pos >cdsStart[gene] and pos < cdsStop[gene]:
						try:
							thisCount[line[2]]+=1
						except KeyError:
							thisCount[line[2]]=1
						c+=1
					individGenes[gene]=pos
				except KeyError: pass
		
		
		for gene in thisCount.keys():
			if thisCount[gene]>readCutoff and cdsStop[gene]-cdsStart[gene]>10:
				rpkm=10**9*thisCount[gene]/(float(c)*(cdsStop[gene]-cdsStart[gene]))
				geneLevel[gene]=rpkm
		
		metageneProfiles=[]	
		for gene in geneLevel.keys():			
			if len(individGenes[gene])>=individReadCutoff:
				thisGeneDist=writePosGraphs(individGenes[gene],0,mRNALen[gene],'plots/individGenes/' + refToCommon[gene] +'-'+ fileName.strip('.map'),gene=gene)
				
				#if geneLevel[gene]>rpkmCutoff:
				if cdsStart[gene]> posxmin*-1 and cdsStop[gene]-cdsStart[gene]>posxmax:	
					metaBuild=thisGeneDist[cdsStart[gene]-20:cdsStart[gene]+posxmax]
					if np.sum(metaBuild)>0:
						metageneProfiles.append(np.array(metaBuild)/np.sum(metaBuild))
					
		#Make metagene profile
		dist=[]
		write=open('plots/riboPosition/' + folder,'w')
		x=range(posxmin,posxmax)
		for pos in range(posxmax-posxmin):
			build=[]
			for gene in metageneProfiles:
				try:
					build.append(gene[pos])
				except KeyError:pass
			dist.append(np.mean(build))
			write.writelines('\n' + str(x[pos]) + '\t' + str(dist[pos]))
		write.close()
		plt.plot(range(posxmin,posxmax),dist,color=(1,0,0))
		plt.xlabel('Position relative to start')
		plt.ylabel('Percentage of reads')
		plt.axis(xmin=posxmin,xmax=posxmax,ymin=0)
		plt.savefig('plots/riboPosition/' + fileName.strip('.map') + '.png')
		plt.savefig('plots/riboPosition/' + fileName.strip('.map') + '.svg')
		plt.clf()
		
		
	except IOError:
		print 'Missing: ' + fileName
	
	return geneLevel
	
	

cdsStart=dict();cdsStop=dict();refToCommon=dict();commToRef=dict();mRNALen=dict()
fullName=dict()	
memProt=dict()
dengueResponse=dict()
allGenes=[]
for line in open('database/geneInfo'):
	if len(line)>4 and line[0] !="#":
		line=line.strip("\n").split("\t")
		allGenes.append(line[0])
		cdsStart[line[0]]=int(line[2]);cdsStop[line[0]]=int(line[3]);mRNALen[line[0]]=int(line[4]);commToRef[line[1]]=line[0];	refToCommon[line[0]]=line[1]

#Load footprints
print 'Printing graphs to the plots/ folder...'
files=os.listdir('mapped/')
for lib in files:
	vars()[lib.strip('.map') + 'abund']=loadReads('mapped',lib)

		
print 'Writing gene expression levels (totalExpression.csv)'
write=open('totalExpression.csv','w')


write.writelines('\t')
for sample in files:
	write.writelines('\t' + sample.strip('.map'))
for gene in allGenes:
	write.writelines('\n' + gene + '\t' + refToCommon[gene])
	for sample in files:
		
		try:
			toWrite='\t' + str(vars()[sample.strip('.map') + 'abund'][gene])
		except KeyError:
			toWrite='\t'
		write.writelines(toWrite)
write.close()

print 'Done! ' + str(len(files)) + ' libraries analyzed.'
