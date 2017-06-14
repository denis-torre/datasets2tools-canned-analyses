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
import sys, os, glob, json, urllib2, xmltodict, gzip, StringIO, requests, tinys3
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
awsFile = 'data/aws.json'

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
########## S1. CREEDS
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

#############################################
########## 2. Submit to G2E
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

#############################################
########## 3. Make Canned Analyses
#############################################

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

	# Filter dataframe
	cannedAnalysisDataframe = mergedDataframe[['dataset_accession', 'tool_name', 'canned_analysis_title', 'canned_analysis_description', 'canned_analysis_url', 'canned_analysis_preview_url', 'metadata']].reset_index(drop=True).dropna()

	# Fix metadata
	cannedAnalysisDataframe['metadata'] = [json.dumps(x) for x in cannedAnalysisDataframe['metadata']]

	# Add date
	cannedAnalysisDataframe['date'] = '2017-05-22'

	# Save file
	cannedAnalysisDataframe.to_csv(outfile, sep='\t', index=False)

#######################################################
#######################################################
########## S2. ARCHS4
#######################################################
#######################################################

#############################################
########## 1. Get JSON
#############################################

def archsJobs():

	# Read link dataframe
	with open('f2-archs4.dir/gse_platform_ids.txt', 'r') as openfile:

		# Read lines
		outfiles = ['f2-archs4.dir/json/'+x.split('_series')[0]+'.json' for x in openfile.readlines()]

		# Get length
		n = len(outfiles)
		i = 0

		# Loop through outfiles
		for outfile in outfiles:

			# Print
			i += 1
			print 'Doing file {i}/{n}...'.format(**locals())

			# Yield
			yield [None, outfile]

@files(archsJobs)

def getArchsJson(infile, outfile):

	try:

		# Get matrix URL
		matrixUrl = 'https://s3.amazonaws.com/mssm-seq-series-platform/'+os.path.basename(outfile)[:-len('.json')]+'_series_matrix.txt.gz'

		# Initialize matrix
		seriesMatrix = SeriesMatrix(matrixUrl)

		# Get expression data
		seriesMatrix.get_expression_data()

		# Get sample metadata
		seriesMatrix.get_sample_metadata()

		# Get JSON
		seriesMatrix.get_clustergrammer_json(outfile)

	except:
		
		# Create file
		os.system('touch {outfile}'.format(**locals()))

#############################################
########## 2. Upload JSON
#############################################

# @follows(getArchsJson)

@transform('f2-archs4.dir/json',
		   suffix(''),
		   add_inputs(awsFile),
		   '_upload.txt')

def uploadJsonFiles(infiles, outfile):

	# Split infiles
	jsonDir, awsFile = infiles

	# Read data
	with open(awsFile, 'r') as openfile:
		aws_dict = json.loads(openfile.read())

	# Define S3 function
	def uploadS3(file, key, bucket):
		conn = tinys3.Connection(aws_dict['awsid'], aws_dict['awskey'], tls=True)
		f = open(file,'rb')
		conn.upload(key, f, bucket)
	
	# Define basename function
	def basename(p):
		temp = p.split("/")
		return temp[len(temp)-1]

	# Upload
	os.chdir(os.path.dirname(jsonDir))
	files = os.listdir("json")
	for f in files:
		uploadS3("json/"+f, basename(f), "mssm-seq-series-json")
		
	# Create file
	os.system('cd ..; touch {outfile}'.format(**locals()))

#############################################
########## 3. Make Canned Analyses
#############################################

@merge(glob.glob('f2-archs4.dir/json/*'),#getArchsJson,
	   'f2-archs4.dir/archs4-canned_analyses.txt')

def makeArchsCannedAnalyses(infiles, outfile):

	# Get analysis dict
	canned_analysis_dict = [{'dataset_accession': os.path.basename(x).split('_')[0], 'tool_name': 'archs4', 'canned_analysis_url':'http://amp.pharm.mssm.edu/datasets2tools/analysis/archs4?q='+os.path.basename(x)[:-len('.json')], 'canned_analysis_title': 'Interactive heatmap visualization of RNA-seq dataset '+os.path.basename(x).split('_')[0], 'canned_analysis_description': 'Highly interactive web-based heatmap visualization of the top 500 most variable genes in the '+os.path.basename(x).split('_')[0]+' RNA-seq dataset, as processed by ARCHS4.', 'canned_analysis_preview_url':'https://github.com/denis-torre/images/blob/master/archs4/'+str(randint(1, 5))+'.png?raw=true', 'metadata': json.dumps({'platform': os.path.basename(x).split('_')[1][:-len('.json')], 'Genes': '500', 'Gene Rank Method': 'Variance', 'Data Normalization Method': 'Z-score'})} for x in infiles]

	# Get dataframe
	canned_analysis_dataframe = pd.DataFrame(canned_analysis_dict)

	# Write
	canned_analysis_dataframe.to_csv(outfile, sep='\t', index=False)

#######################################################
#######################################################
########## S3. GeneMANIA
#######################################################
#######################################################

#############################################
########## 1. Process analyses
#############################################

@follows(mkdir('f3-genemania.dir'))

@transform(getCreedsData,
		   regex(r'.*/(.*).csv'),
		   r'f3-genemania.dir/\1-genemania_links.txt')

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


##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')
