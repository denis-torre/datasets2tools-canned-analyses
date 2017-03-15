#################################################################
#################################################################
############### Canned Analyses ################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
from ruffus import *
import sys, os, glob, json, urllib2, xmltodict
import pandas as pd
import rpy2.robjects as robjects
import pandas.rpy.common as com
import numpy as np
from lxml import etree

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
sys.path.append('/Users/denis/Documents/Projects/creeds/creeds_orm')
import Support as S
import PipelineDatasets2toolsCannedAnalyses as P
import orm

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
metadataFileNames = ['single_gene_perturbations-v1.0.csv', 'disease_signatures-v1.0.csv', 'single_drug_perturbations-v1.0.csv', 'single_drug_perturbations-DM.csv', 'single_gene_perturbations-p1.0.csv', 'disease_signatures-p1.0.csv', 'single_drug_perturbations-p1.0.csv']
dbFile = 'g2e/db.json'

##### 2. R Connection #####
rSource = 'pipeline/scripts/pipeline-datasets2tools-canned-analyses.R'
r = robjects.r
r.source(rSource)

#################################################################
#################################################################
############### 1. CREEDS #######################################
#################################################################
#################################################################

#######################################################
#######################################################
########## S1. Download Data
#######################################################
#######################################################

#############################################
########## 1. Download metadata
#############################################

@follows(mkdir('f1-creeds.dir'))

def downloadJobs():

	# Loop through metadata
	for metadataFileName in metadataFileNames:

		# Get dirname
		dirName = os.path.join('f1-creeds.dir', metadataFileName.replace('.csv', ''))

		# Create
		if not os.path.exists(dirName):
			os.makedirs(dirName)

		# Get outfile
		outfile = os.path.join(dirName, metadataFileName)

		# Yield
		yield [None, outfile]

@files(downloadJobs)

def getCreedsData(infile, outfile):

	# Get filename
	fileName = os.path.basename(outfile)

	# Get link
	metadataFileLink = 'http://amp.pharm.mssm.edu/CREEDS/download/%(fileName)s' % locals()

	# Get outfile directory
	outDir = os.path.dirname(outfile)

	# Prepare statement
	os.system('wget -P %(outDir)s %(metadataFileLink)s' % locals())

#######################################################
#######################################################
########## S2. GEO2Enrichr
#######################################################
#######################################################

#############################################
########## 1. Submit to G2E
#############################################

@transform(getCreedsData,
		   regex(r'(.*)/.*/(.*).csv'),
		   r'\1/\2/\2-links.txt')

def getCreedsLinks(infile, outfile):

	# Report
	print 'Doing', infile, '...'

	# Get dataframe
	creedsDataframe = pd.read_csv(infile, index_col='id').iloc[:10,]

	# Initialize dict
	linkDict = {x:{} for x in creedsDataframe.index}

	# Loop through IDs
	for creedsId in creedsDataframe.index:

		print creedsId

		# Get DBSignature
		sig = orm.DBSignature(creedsId)

		# Initialize vectors
		sig.init_cs_vectors(cutoff=2000)

		# Get PAEA
		try:
			linkDict[creedsId]['paea'] = sig.post_to_paea()
		except:
			linkDict[creedsId]['paea'] = None

		# Get L1000
		linkDict[creedsId]['l1000cds2'] = sig.post_to_cds2()


		# Get enrichr
		enrichrDict = sig.post_to_enrichr()
		linkDict[creedsId]['enrichr-up'] = enrichrDict['up']
		linkDict[creedsId]['enrichr-down'] = enrichrDict['down']

		# Convert to dataframe
		linkDataframe = pd.DataFrame(linkDict).T

	# Reset index
	linkDataframe['creeds_id'] = linkDataframe.index

	# Melt
	linkDataframeMelt = pd.melt(linkDataframe, id_vars='creeds_id', var_name='analysis', value_name='link')

	# Save file
	linkDataframeMelt.to_csv(outfile, sep='\t', index=False)

#######################################################
#######################################################
########## S2. Canned Analyses
#######################################################
#######################################################

@transform(getCreedsLinks,
		   regex(r'(.*)/.*/(.*)-links.txt'),
		   add_inputs(r'\1/\2/\2.csv'),
		   r'\1/\2/\2-canned_analyses.txt')

def makeCreedsCannedAnalyses(infiles, outfile):

	# Split infiles
	linkFile, metadataFile = infiles

	# Read dataframes
	linkDataframe = pd.read_table(linkFile, index_col='creeds_id')
	metadataDataframe = pd.read_csv(metadataFile, index_col='id')

	# Process dataframes
	linkDataframe['tool'] = [x.split('-')[0] if '-' in x else x for x in linkDataframe['analysis']]
	linkDataframe['geneset'] = [x.split('-')[1]+'regulated' if '-' in x else 'combined' for x in linkDataframe['analysis']]
	linkDataframe['top_genes'] = [str(500) if 'enrichr' in x else np.nan for x in linkDataframe['analysis']]
	metadataDataframe['ctrl_ids'] = [x.replace('|', ', ') for x in metadataDataframe['ctrl_ids']]
	metadataDataframe['pert_ids'] = [x.replace('|', ', ') for x in metadataDataframe['pert_ids']]

	# Merge metadata dataframe
	mergedDataframe = pd.merge(metadataDataframe, linkDataframe, left_index=True, right_index=True, how='inner')

	# Add metadata column
	mergedMetadataDataframe = mergedDataframe.drop(['geo_id', 'platform', 'version', 'link', 'analysis', 'tool'], axis=1)
	mergedDataframe['metadata'] = [json.dumps(dict(rowData.dropna())) for index, rowData in mergedMetadataDataframe.iterrows()]

	# Filter dataframe
	cannedAnalysisDataframe = mergedDataframe[['geo_id', 'tool', 'link', 'metadata']].reset_index(drop=True)

	# Save file
	cannedAnalysisDataframe.to_csv(outfile, sep='\t', index_label='index')

#################################################################
#################################################################
############### 2. ARCHS4 #######################################
#################################################################
#################################################################

#######################################################
#######################################################
########## S1. ARCHS4
#######################################################
#######################################################

#############################################
########## 1. Read links
#############################################

@follows(mkdir('f2-archs.dir'))

@files(None,
	   'f2-archs.dir/archs-links.txt')

def getArchsLinks(infile, outfile):

	# Read XML
	xmlString = urllib2.urlopen('https://s3.amazonaws.com/mssm-seq-series').read()

	# Convert to dict
	xmlDict = xmltodict.parse(xmlString)

	# Get GEO IDs
	geoIds = [str(x['Key'].split('_')[0]) for x in xmlDict['ListBucketResult']['Contents']]

	# Get data
	linkDict = {x: json.loads(urllib2.urlopen('http://amp.pharm.mssm.edu/awsscheduler/getGSEmatrix.php?gse=%(x)s' % locals()).read()) for x in geoIds}

	# Convert to dataframe
	linkDataframe = pd.DataFrame(linkDict).T

	# Write file
	linkDataframe.to_csv(outfile, sep='\t', index=False)

#################################################################
#################################################################
############### 3. GEO ##########################################
#################################################################
#################################################################

#######################################################
#######################################################
########## S1. Get data
#######################################################
#######################################################

#############################################
########## 1. Get GSE
#############################################

def gseJobs():

	for organism in ['human', 'mouse']:

		outfile = 'f3-geo.dir/%(organism)s/%(organism)s-geo_files.txt' % locals()

		dirName = os.path.dirname(outfile)

		if not os.path.exists(dirName):
			os.makedirs(dirName)

		yield [None, outfile]

@follows(mkdir('f3-geo.dir'))

@files(gseJobs)

def getGseFiles(infile, outfile):

	# Get organism
	organism = os.path.basename(outfile).split('-')[0]

	# Get eSearch URL
	eSearchURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=(Expression%20profiling%20by%20array%5BDataSet%20Type%5D)%20AND%20' + organism + '%5BOrganism%5D&retMax=100000'

	# Read data
	idList = xmltodict.parse(urllib2.urlopen(eSearchURL).read())['eSearchResult']['IdList']['Id']

	# Initialize result
	resultList = []

	# Split in subsets
	for idListSubset in np.array_split(idList, len(idList)/50):

	    # Create comma-separated string
	    idString = ','.join(idListSubset)

	    # Get eSummary URL
	    eSummaryURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=%(idString)s' % locals()

	    # Get root
	    root = etree.fromstring(urllib2.urlopen(eSummaryURL).read())
	    
	    # Get data
	    elemData = [[y.text for y in x.getchildren() if y.get('Name') in ['Accession', 'GDS', 'FTPLink']] for x in root]
	    
	    # Append
	    resultList += elemData

	# Convert to dataframe
	resultDataframe = pd.DataFrame(resultList, columns=['gds', 'id', 'ftp_link'])

	# Fix FTP link
	resultDataframe['ftp_link'] = [os.path.join(ftp_link, 'soft', gds+'.soft.gz') for gds, ftp_link in resultDataframe[['gds', 'ftp_link']].as_matrix()]

	# Save file
	resultDataframe.to_csv(outfile, sep='\t', index=False)

##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')
