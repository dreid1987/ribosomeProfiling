import os



def which(program): #From http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
	def is_exe(fpath):
		return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

	fpath, fname = os.path.split(program)
	if fpath:
		if is_exe(program):
			return program
	else:
 		for path in os.environ["PATH"].split(os.pathsep):
			path = path.strip('"')
			exe_file = os.path.join(path, program)
			if is_exe(exe_file):
				return exe_file

	return ""
def addNextFolder(currentFolder):
	folders=os.listdir(currentFolder)
	return currentFolder + folders[0] +'/'



print 'Checking dependencies...'
print 'You may be asked to enter your OS password if dependencies need to be installed'

deps=dict()
deps['pip']='python-pip'
deps['samtools']='samtools'
deps['bowtie']='bowtie'
deps['cufflinks']='cufflinks'
deps['tophat']='tophat'

for dep in deps.keys():
	if len(which(dep))==0:
		os.system('sudo apt-get install ' + deps[dep])

if len(which('cutadapt'))==0:
	os.system('sudo pip install cutadapt')
try:
	import matplotlib.pyplot as plt
except ImportError:
	os.system('sudo pip install matplotlib')
try:
	import Bio
except ImportError:
	os.system('sudo apt-get install python-biopython')


print 'What is the name of your reference genome file?'
print 'Before running this program, download a Refseq genome and extract it into a folder called \'database\' that is inside your currect directory'
print 'You may download these from http://ccb.jhu.edu/software/tophat/igenomes.shtml.'
print 'Example input: Homo_sapiens'
refFolder = raw_input('> ')



#Find database
databaseFolder='database/' + refFolder + '/'
databaseFolder=addNextFolder(databaseFolder)
databaseFolder=addNextFolder(databaseFolder)
	
gtfOK = os.path.isfile(databaseFolder + 'Annotation/Genes/genes.gtf') 
if gtfOK:
	gtfLoc=databaseFolder + 'Annotation/Genes/genes.gtf'
else:	
	while not gtfOK:
		gtfLoc=raw_input('GTF file not found. Specify GTF file: ')
		gtfOK = os.path.isfile(gtfLoc) 
chrOK = os.path.isfile(databaseFolder + "Sequence/Chromosomes/1.fa") 
if chrOK:
	chrLoc=databaseFolder + "Sequence/Chromosomes/"
else:	
	while not chrOK:
		chrLoc=raw_input('Genome sequence files not found. Specify folder with .fa files: ')	
		chrOK = os.path.isfile(chrLoc + '1.fa') 

if not os.path.isdir('trim'):
	os.mkdir('trim')
if not os.path.isdir('mapped'):
	os.mkdir('mapped')



readFiles=os.listdir('reads/')
print 'What is the 3\' adapter to trim?'
print "Enter 'N' to use the default NEB or Illumina adapters, or enter your own sequence"
toTrim = raw_input('> ')
if toTrim=='N':
	toTrim='AAGATCGGAAGAGCACACGTCT'

print '\n Trimming adapters from reads using Cutadapt'	
for readFile in readFiles:
	os.system('cutadapt -a ' + toTrim +' --discard-untrimmed -m 20 reads/' + readFile + ' -o trim/' + readFile)

print '\n Building a transcriptome using Tophat Cufflinks'	
tophatCommand= 'tophat -G ' + gtfLoc + ' -T -p 3 --bowtie1 ' + databaseFolder + 'Sequence/BowtieIndex/genome '
readFilePos=1
for readFile in readFiles:
	tophatCommand = tophatCommand + 'trim/' + readFile
	if readFilePos<len(readFiles):
		tophatCommand = tophatCommand +','
os.system(tophatCommand)
		
os.system('cufflinks -m 35 -s 5 -G ' + databaseFolder + 'Annotation/Genes/genes.gtf tophat_out/accepted_hits.bam')



os.system('python chooseTranscripts.py ' + gtfLoc + ' ' + chrLoc)
os.system('bowtie-build database/chosenTranscripts database/chosenTranscripts')
for readFile in readFiles:
	os.system('bowtie -p 3 --norc -v 0 database/chosenTranscripts trim/' + readFile + ' mapped/' + readFile.split('.')[0] + '.map')



if not os.path.isdir('plots'):
	os.system('mkdir plots')
if not os.path.isdir('plots/riboPosition'):
	os.system('mkdir plots/riboPosition')
if not os.path.isdir('plots/individGenes'):
	os.system('mkdir plots/individGenes')

os.system('python diff.py')
