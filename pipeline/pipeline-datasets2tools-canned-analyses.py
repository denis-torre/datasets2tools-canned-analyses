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
import sys, os, glob, json, urllib2, xmltodict, gzip, StringIO, requests
import pandas as pd
import rpy2.robjects as robjects
# import pandas.rpy.common as com
import numpy as np
from lxml import etree
from collections import Counter
from random import randint

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
sys.path.append('/Users/denis/Documents/Projects/creeds/creeds_orm')
import Support as S
import PipelineDatasets2toolsCannedAnalyses as P
import orm
from SeriesMatrix import *

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
metadataFileNames = ['single_gene_perturbations-v1.0.csv', 'disease_signatures-v1.0.csv', 'single_drug_perturbations-v1.0.csv']#, 'single_drug_perturbations-DM.csv', 'single_gene_perturbations-p1.0.csv', 'disease_signatures-p1.0.csv', 'single_drug_perturbations-p1.0.csv']
dbFile = 'data/db.json'
screenshotFile = 'data/screenshots/screenshots.txt'

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
	creedsDataframe = pd.read_csv(infile, index_col='id')

	# Initialize dict
	linkDict = {x:{} for x in creedsDataframe.index}

	i = 0
	n = len(creedsDataframe.index)

	# Loop through IDs
	for creedsId in creedsDataframe.index:

		i += 1
		print 'Doing %(creedsId)s, %(i)s/%(n)s...' % locals()

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
		try:
			linkDict[creedsId]['l1000cds2-mimic'] = sig.post_to_cds2(aggravate=True)
		except:
			linkDict[creedsId]['l1000cds2-mimic'] = None
		try:
			linkDict[creedsId]['l1000cds2-reverse'] = sig.post_to_cds2(aggravate=False)
		except:
			linkDict[creedsId]['l1000cds2-reverse'] = None

		# Get enrichr
		try:
			enrichrDict = sig.post_to_enrichr()
			linkDict[creedsId]['enrichr-up'] = enrichrDict['up']
			linkDict[creedsId]['enrichr-down'] = enrichrDict['down']
		except:
			linkDict[creedsId]['enrichr-up'] = None
			linkDict[creedsId]['enrichr-down'] = None

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
		   add_inputs(r'\1/\2/\2.csv',
		   			  screenshotFile),
		   r'\1/\2/\2-canned_analyses.txt')

def makeCreedsCannedAnalyses(infiles, outfile):

	# Split infiles
	linkFile, metadataFile, screenshotFile = infiles

	# Read dataframes
	linkDataframe = pd.read_table(linkFile, index_col='creeds_id').rename(columns={'link': 'canned_analysis_url'})
	metadataDataframe = pd.read_csv(metadataFile).rename(columns={'geo_id': 'dataset_accession', 'id': 'creeds_id'}).set_index('creeds_id', drop=False)
	toolScreenshotDataframe = pd.read_table(screenshotFile, index_col='screenshot_label')

	# Process dataframes
	linkDataframe['tool_name'] = [x.split('-')[0] if '-' in x else x for x in linkDataframe['analysis']]
	linkDataframe['geneset'] = [x.split('-')[1]+'regulated' if 'enrichr-' in x else np.nan for x in linkDataframe['analysis']]
	linkDataframe['direction'] = [x.split('-')[1] if 'l1000cds2-' in x else np.nan for x in linkDataframe['analysis']]
	linkDataframe['top_genes'] = [str(500) if 'enrichr' in x else np.nan for x in linkDataframe['analysis']]
	linkDataframe['canned_analysis_preview_url'] = [toolScreenshotDataframe.loc[x].ix[randint(0, 2), 'screenshot_url'] for x in linkDataframe['analysis']]
	metadataDataframe['ctrl_ids'] = [x.replace('|', ', ') for x in metadataDataframe['ctrl_ids']]
	metadataDataframe['pert_ids'] = [x.replace('|', ', ') for x in metadataDataframe['pert_ids']]
	metadataDataframe['cell_type'] = [x.replace('"', '') if type(x) == str else x for x in metadataDataframe['cell_type']]

	# Merge metadata dataframe
	mergedDataframe = pd.merge(metadataDataframe, linkDataframe, left_index=True, right_index=True, how='inner').drop(['platform', 'version', 'analysis'], axis=1)

	# Add metadata column
	def removeNonAscii(s): return "".join([x if ord(x) < 128 else '?' for x in s])
	mergedDataframe['metadata'] = [{key:value if type(value) == float else removeNonAscii(value).replace('"', '') for key, value in rowData.dropna().iteritems()} for index, rowData in mergedDataframe.drop(['canned_analysis_preview_url', 'dataset_accession', 'canned_analysis_url', 'tool_name'], axis=1).iterrows()]

	# Add description and title
	mergedDataframe['canned_analysis_title'] = [P.getAnalysisTitle(tool_name, metadata) for tool_name, metadata in mergedDataframe[['tool_name', 'metadata']].as_matrix()]
	mergedDataframe['canned_analysis_description'] = [P.getAnalysisDescription(tool_name, metadata) for tool_name, metadata in mergedDataframe[['tool_name', 'metadata']].as_matrix()]

	# JSON metadata
	mergedDataframe['metadata'] = [json.dumps(x) for x in mergedDataframe['metadata']]

	# Filter dataframe
	cannedAnalysisDataframe = mergedDataframe[['dataset_accession', 'tool_name', 'canned_analysis_title', 'canned_analysis_description', 'canned_analysis_url', 'canned_analysis_preview_url', 'metadata']].reset_index(drop=True)

	# Save file
	cannedAnalysisDataframe.to_csv(outfile, sep='\t', index=False)

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
	linkList = []

	# Split in subsets
	for idListSubset in np.array_split(idList, len(idList)/100):

		# Create comma-separated string
		idString = ','.join(idListSubset)

		# Get eSummary URL
		eSummaryURL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=%(idString)s' % locals()

		# Get root
		root = etree.fromstring(urllib2.urlopen(eSummaryURL).read())
		
		# Get data
		elemData = [[y.text for y in x.getchildren() if y.get('Name') in ['Accession', 'FTPLink', 'GPL']] for x in root]
		
		# Append
		linkList += elemData
		
	# Convert to dataframe
	linkDataframe = pd.DataFrame(linkList, columns=['geo_accession', 'gpl', 'ftp_link'])

	# Get result list
	resultList = [P.getMatrixLink(geo_accession, gpl, ftp_link) for geo_accession, gpl, ftp_link in linkDataframe[['geo_accession', 'gpl', 'ftp_link']].as_matrix() if gpl != None]
	resultList = [item for sublist in resultList for item in sublist]

	# Get link dataframe
	resultDataframe = pd.DataFrame(resultList, columns=['geo_accession', 'platform' , 'matrix_link'])

	# Write
	resultDataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Split data
#############################################

@subdivide(getGseFiles,
		   regex(r'(.*)/.*-geo_files.txt'),
		   r'\1/*/*-files.txt',
		   r'\1')

def splitGeoFiles(infile, outfiles, outfileRoot):

	# Read data
	geoDataframe = pd.read_table(infile).set_index('platform', drop=False)

	# Get counter
	platformCounts = Counter(geoDataframe.index)

	# Select platforms
	platforms = [platform for platform, count in platformCounts.iteritems() if count > 500]

	# Loop
	for platform in platforms:

		# Get subset
		geoDataframeSubset = geoDataframe.loc[platform]

		# get outfile
		outfile = "{outfileRoot}/{platform}/{platform}-files.txt".format(**locals())

		# Get directory name
		outDir = os.path.dirname(outfile)

		# Create directory, if not exists
		if not os.path.exists(outDir):
			os.makedirs(outDir)

		# Write
		geoDataframeSubset.to_csv(outfile, sep='\t', index=False)

#############################################
########## 3. Get links
#############################################

@transform(splitGeoFiles,
		   regex(r'(.*)/(.*)/(.*)/.*-files.txt'),
		   add_inputs(r'\1/platform_annotations/\3-annotation.txt'),
		   r'\1/\2/\3/\3-links.txt')

def getClustergramLinks(infiles, outfile):

	# Split infiles
	geoFile, platformFile = infiles

	# Read data
	geoFileDataframe = pd.read_table(geoFile)

	# Initialize result list
	resultList = []

	i = 0
	n = len(geoFileDataframe.index)
	organism = geoFile.split('/')[1]

	# Loop through dataframe
	for geo_accession, platform, matrix_link in geoFileDataframe.as_matrix():

		i += 1

		print 'Doing sample %(i)s/%(n)s for %(organism)s platform %(platform)s...' % locals()

		try:

			# Get response
			response = urllib2.urlopen(matrix_link)

			# Read file
			geoData = [x.strip() for x in gzip.GzipFile(fileobj=StringIO.StringIO(response.read())).readlines()]

			# If soft
			if 'soft' in matrix_link:
				clustergramDataframe = P.prepareSoftFile(geoData)
			elif 'series' in matrix_link:
				clustergramDataframe = P.prepareSeriesMatrix(geoData, platformFile)
			else:
				raise ValueError('Unknown type: ' + matrix_link)

			# Zscore
			clustergramDataframe = ((clustergramDataframe.T - clustergramDataframe.T.mean())/clustergramDataframe.T.std()).T

			# Get filename
			tempfile = os.path.basename(matrix_link)+'.txt'

			# Save file
			clustergramDataframe.to_csv(tempfile, sep='\t', index_label='Gene Symbol')

			# Get link
			upload_url = 'http://amp.pharm.mssm.edu/clustergrammer/matrix_upload/'

			# Post request
			r = requests.post(upload_url, files={'file': open(tempfile, 'rb')})

			# Remove file
			os.unlink(tempfile)

			# Get link string
			link = os.path.dirname(r.text)

		except:

			# Set error
			link = None

		# Add to dict
		resultList.append(link)

	# Convert to dataframe
	geoFileDataframe['link'] = resultList

	# Save
	geoFileDataframe.to_csv(outfile, sep='\t', index_label='geo_accession')

#################################################################
#################################################################
############### 4. LINCS ########################################
#################################################################
#################################################################

#######################################################
#######################################################
########## S1. Process Data
#######################################################
#######################################################

#############################################
########## 1. Process analyses
#############################################

@follows(mkdir('f4-lincs.dir'))

@transform('f4-lincs.dir/lincs_canned_analyses.csv',
		   suffix('.csv'),
		   '.txt')

def processLincsAnalyses(infile, outfile):
	
	# Read dataframe
	analysisDataframe = pd.read_csv(infile).set_index('dataset_accessions_list', drop=False)

	# Get subset
	cannedAnalysisDataframe = analysisDataframe.rename(columns={'title': 'canned_analysis_title', 'dataset_accessions_list': 'dataset_accession', 'ca_image_url': 'canned_analysis_preview_url'}).loc['LDG-1285', ['dataset_accession', 'tool_name', 'canned_analysis_title', 'canned_analysis_description', 'canned_analysis_url', 'canned_analysis_preview_url']]

	# Get metadata
	cannedAnalysisDataframe['metadata'] = [json.dumps({'cell_line': 'A375, A549, MCF7, NPC, PC3', 'perturbation_type': 'epigenetic compound'}),
										   json.dumps({'cell_line': 'MCF7, NPC', 'perturbation_type': 'kinase inhibitors'})]

	# Save
	cannedAnalysisDataframe.to_csv(outfile, sep='\t', index=False)

#################################################################
#################################################################
############### 5. GeneMANIA ####################################
#################################################################
#################################################################

#######################################################
#######################################################
########## S1. Create Analyses
#######################################################
#######################################################

#############################################
########## 1. Process analyses
#############################################

@follows(mkdir('f5-genemania.dir'))

@transform(getCreedsData,
		   regex(r'.*/(.*).csv'),
		   r'f5-genemania.dir/\1-genemania_links.txt')

def getGenemaniaLinks(infile, outfile):

	# Print
	print 'Doing '+infile+'...'

	# Get dataframe
	creedsDataframe = pd.read_csv(infile, index_col='id')

	# Initialize dict
	linkDict = {x:{} for x in creedsDataframe.index}

	# Loop through dataframe
	for creedsId, rowData in creedsDataframe.iterrows():
		
		# Get DBSignature
		sig = orm.DBSignature(creedsId)

		# Initialize vectors
		sig.init_cs_vectors(cutoff=2000)
		
		# Link dict
		try:
			linkDict[creedsId] = sig.get_genemania_url(rowData['organism'])
		except:
			linkDict[creedsId] = {'up': None, 'down': None}

	# Convert to dataframe
	linkDataframe = pd.DataFrame(linkDict).T

	# Reset index
	linkDataframe['creeds_id'] = linkDataframe.index

	# Melt
	linkDataframeMelt = pd.melt(linkDataframe, id_vars='creeds_id', var_name='analysis', value_name='link')

	# Save
	linkDataframeMelt.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Get canned analyses
#############################################

@transform(getGenemaniaLinks,
		   regex(r'(.*)/(.*)-genemania_links.txt'),
		   add_inputs(r'f1-creeds.dir/\2/\2.csv'),
		   r'\1/\2-genemania_canned_analyses.txt')

def makeGenemaniaCannedAnalyses(infiles, outfile):

	# Split infiles
	linkFile, metadataFile = infiles

	# Read dataframes
	linkDataframe = pd.read_table(linkFile, index_col='creeds_id').rename(columns={'link': 'canned_analysis_url'})
	metadataDataframe = pd.read_csv(metadataFile).rename(columns={'geo_id': 'dataset_accession', 'id': 'creeds_id'}).set_index('creeds_id', drop=False)
	toolScreenshotDataframe = pd.read_table(screenshotFile, index_col='screenshot_label')

	# Process dataframes
	linkDataframe['tool_name'] = 'genemania'
	linkDataframe['geneset'] = [x+'regulated' for x in linkDataframe['analysis']]
	linkDataframe['top_genes'] = str(50)
	linkDataframe['canned_analysis_preview_url'] = 'https://raw.githubusercontent.com/denis-torre/images/master/genemania/'+str(randint(1, 5))+'.png'
	metadataDataframe['ctrl_ids'] = [x.replace('|', ', ') for x in metadataDataframe['ctrl_ids']]
	metadataDataframe['pert_ids'] = [x.replace('|', ', ') for x in metadataDataframe['pert_ids']]
	metadataDataframe['cell_type'] = [x.replace('"', '') if type(x) == str else x for x in metadataDataframe['cell_type']]

	# Merge metadata dataframe
	mergedDataframe = pd.merge(metadataDataframe, linkDataframe, left_index=True, right_index=True, how='inner').drop(['platform', 'version', 'analysis'], axis=1)

	# Add metadata column
	def removeNonAscii(s): return "".join([x if ord(x) < 128 else '?' for x in s])
	mergedDataframe['metadata'] = [{key:value if type(value) == float else removeNonAscii(value).replace('"', '') for key, value in rowData.dropna().iteritems()} for index, rowData in mergedDataframe.drop(['canned_analysis_preview_url', 'dataset_accession', 'canned_analysis_url', 'tool_name'], axis=1).iterrows()]

	# Add description and title
	mergedDataframe['canned_analysis_title'] = [P.getAnalysisTitle(tool_name, metadata) for tool_name, metadata in mergedDataframe[['tool_name', 'metadata']].as_matrix()]
	mergedDataframe['canned_analysis_description'] = [P.getAnalysisDescription(tool_name, metadata) for tool_name, metadata in mergedDataframe[['tool_name', 'metadata']].as_matrix()]

	# JSON metadata
	mergedDataframe['metadata'] = [json.dumps(x) for x in mergedDataframe['metadata']]

	# Filter dataframe
	cannedAnalysisDataframe = mergedDataframe[['dataset_accession', 'tool_name', 'canned_analysis_title', 'canned_analysis_description', 'canned_analysis_url', 'canned_analysis_preview_url', 'metadata']].reset_index(drop=True)

	# Save file
	cannedAnalysisDataframe.to_csv(outfile, sep='\t', index=False)

#################################################################
#################################################################
############### 6. ARCHS4 #######################################
#################################################################
#################################################################

#######################################################
#######################################################
########## S1. Create Analyses
#######################################################
#######################################################

#############################################
########## 1. Get links
#############################################

@follows(mkdir('f6-archs4.dir'))

@originate('f6-archs4.dir/archs4-matrix_links.txt')

def getMatrixLinks(outfile):

	# Get XML
	xmlDict = xmltodict.parse(urllib2.urlopen('https://s3.amazonaws.com/mssm-seq-series').read())

	# Get accessions
	geoAccessions = [x['Key'].split('_')[0] for x in xmlDict['ListBucketResult']['Contents']]

	# Get dataframe
	linkDataframe = pd.DataFrame({'dataset_accession': geoAccessions,
	                              'series_matrix_url': ['https://s3.amazonaws.com/mssm-seq-series/{x}_series_matrix.txt.gz'.format(**locals()) for x in geoAccessions]})

	# Save
	linkDataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Get JSON
#############################################

# def vizJobs():

	# Read link dataframe
	# linkDataframe = pd.read_table('f6-archs4.dir/archs4-matrix_links.txt')


@follows(mkdir('f6-archs4.dir/json'))

@subdivide(getMatrixLinks,
		   regex(r'(.*)/archs4-matrix_links.txt'),
		   r'\1/json/*-viz.json',
		   r'\1/json')

def getVizJson(infile, outfiles, outfileRoot):

	# Read link dataframe
	linkDataframe = pd.read_table(infile)

	# Loop through dataframe
	for dataset_accession, series_matrix_url in linkDataframe.iloc[:5].as_matrix():

		# Initialize series matrix
		seriesMatrix = SeriesMatrix(series_matrix_url)

		# Get expression data
		seriesMatrix.get_expression_data()

		# Get sample metadata
		seriesMatrix.get_sample_metadata()

		# Get outfile
		outfile = '{outfileRoot}/{dataset_accession}.json'.format(**locals())

		# Write JSON
		seriesMatrix.get_clustergrammer_json(outfile)


##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')
