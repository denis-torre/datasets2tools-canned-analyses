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
from sqlalchemy import create_engine

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
import Support as S
import PipelineDatasets2toolsCannedAnalyses as P

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

#############################################
########## 2. Merge and filter data
#############################################

@transform(getCreedsData,
		   regex(r'(.*)/.*/(.*).csv'),
		   r'\1/\2/\2_filtered.txt')

def filterCreedsData(infile, outfile):

	# Read dataframe
	creedsDataframe = pd.read_csv(infile, usecols=['id', 'ctrl_ids', 'pert_ids', 'organism', 'geo_id', 'platform'], index_col='id')

	# Get IDs of rows to drop (more than one organism and/or platform)
	rowsToDrop = [creedsId for creedsId, rowData in creedsDataframe.iterrows() if '[' in rowData['organism']+rowData['platform']]

	# Drop
	creedsDataframe.drop(rowsToDrop, inplace=True)

	# Save file
	creedsDataframe.to_csv(outfile, sep='\t', index_label='creeds_id')

#############################################
########## 3. Merge and filter metadata
#############################################

@transform(getCreedsData,
		   regex(r'(.*)/.*/(.*).csv'),
		   add_inputs(r'\1/\2/\2_filtered.txt'),
		   r'\1/\2/\2-metadata.txt')

def processCreedsMetadata(infiles, outfile):

	# Split infiles
	dataFile, filteredDataFile = infiles

	# Read CREEDS dataframe
	metadataDataframe = pd.read_csv(dataFile).set_index('id', drop=False)

	# Read filtered IDs
	filteredCreedsIds = pd.read_table(filteredDataFile)['creeds_id'].tolist()

	# Filter
	metadataDataframe = metadataDataframe.drop(['geo_id', 'platform', 'ctrl_ids', 'pert_ids', 'version'], axis=1).loc[filteredCreedsIds,:]

	# Melt
	metadataDataframeMelt = pd.melt(metadataDataframe, id_vars='id').dropna()

	# Save data
	metadataDataframeMelt.to_csv(outfile, sep='\t', index=False)


#######################################################
#######################################################
########## S2. GEO2Enrichr
#######################################################
#######################################################

#############################################
########## 1. Submit to G2E
#############################################

@transform('f1-creeds.dir/single_drug_perturbations-v1.0/single_drug_perturbations-v1.0_filtered.txt',
		   regex(r'(.*)/.*/(.*)_filtered.txt'),
		   r'\1/\2/\2-extraction_ids.txt')

def submitGeo2Enrichr(infile, outfile):

	# Get dataframe
	creedsDataframe = pd.read_table(infile, index_col='creeds_id').iloc[:10,]

	# Replace with lists
	creedsDataframe['pert_ids'] = [x.split('|') for x in creedsDataframe['pert_ids']]
	creedsDataframe['ctrl_ids'] = [x.split('|') for x in creedsDataframe['ctrl_ids']]

	# Initialize result dict
	resultDict = {}

	# Initialize I
	i = 0
	tot = len(creedsDataframe.index)

	# Loop through rows
	for creedsId, rowData in creedsDataframe.iterrows():

		# Counter
		i += 1
		print 'Doing dataset %(i)s/%(tot)s...' % locals()

		# Submit to GEO2Enrichr
		extractionId = S.submitG2E(dataset = rowData['geo_id'],
								   platform = rowData['platform'],
								   A_cols = rowData['pert_ids'],
								   B_cols = rowData['ctrl_ids'])

		# Add to dict
		resultDict[creedsId] = extractionId

	# Convert to dataframe
	resultDataframe = pd.DataFrame.from_dict(resultDict, orient='index').rename(columns={0: 'extraction_id'})

	# Save file
	resultDataframe.to_csv(outfile, sep='\t', index_label='creeds_id')

#############################################
########## 2. Get links
#############################################

@transform(submitGeo2Enrichr,
		   regex(r'(.*)/.*/(.*)-extraction_ids.txt'),
		   add_inputs(dbFile),
		   r'\1/\2/\2-links.txt')

def getCreedsLinks(infiles, outfile):

	# Split infiles
	extractionFile, dbFile = infiles

	# Get extraction IDs
	extractionIds = pd.read_table(extractionFile)['extraction_id'].tolist()

	# Get string
	extractionIdsString = "('" + "', '".join(extractionIds) + "')"

	# Create DB engine
	with open(dbFile) as openfile:
	    engine = create_engine(os.path.join(json.load(openfile)['phpmyadmin'], 'euclid'))

	# Create SQL command
	sqlCommand = '''SELECT gs.extraction_id AS extraction_id, gl.direction, ta.name, tal.link
	                    FROM gene_signature gs
	                        LEFT JOIN gene_list gl
	                        ON gs.id = gl.gene_signature_fk
	                            LEFT JOIN target_app_link tal
	                            ON gl.id = tal.gene_list_fk
	                            	LEFT JOIN target_app ta
	                            	ON ta.id = tal.target_app_fk
	                    WHERE gs.extraction_id IN %(extractionIdsString)s ''' % locals()

	# Read links
	linkDataframe = pd.read_sql_query(sqlCommand, engine)

	# Save file
	linkDataframe.to_csv(outfile, sep='\t', index=False)

#######################################################
#######################################################
########## S3. Canned Analyses
#######################################################
#######################################################

@transform(getCreedsLinks,
		   regex(r'(.*)/.*/(.*)-links.txt'),
		   add_inputs(r'\1/\2/\2-extraction_ids.txt',
		   			  r'\1/\2/\2.csv'),
		   r'\1/\2/\2-canned_analyses.txt')

def makeCreedsCannedAnalyses(infiles, outfile):

	# Split infiles
	linkFile, idFile, metadataFile = infiles

	# Read dataframes
	linkDataframe = pd.read_table(linkFile)
	idDataframe = pd.read_table(idFile)
	metadataDataframe = pd.read_csv(metadataFile, index_col='id')

	# Process dataframes
	linkDataframe['geneset'] = [{-1: 'underexpressed', 0: 'combined', 1: 'overexpressed'}[x] for x in linkDataframe['direction']]
	metadataDataframeFiltered = metadataDataframe.drop(['geo_id', 'platform', 'ctrl_ids', 'pert_ids', 'version'], axis=1)

	# Create metadata dict
	metadataDict = {creeds_id: {variable: value for variable, value in metadataDataframeFiltered.loc[creeds_id].iteritems() if value == value and value != None} for creeds_id in metadataDataframeFiltered.index}

	# Merge dataframes
	mergedDataframe = idDataframe.merge(linkDataframe, on='extraction_id', how='left')

	# Loop through rows
	for index, rowData in mergedDataframe.iterrows():

		# Get metadata dict
		analysisMetadataDict = metadataDict[rowData['creeds_id']]

		# Add direction
		analysisMetadataDict['geneset'] = rowData['geneset']

		# Add diff exp
		analysisMetadataDict['diff_exp_method'] = 'chdir'
		analysisMetadataDict['genes'] = 500

		# Convert to JSON and add
		mergedDataframe.loc[index, 'metadata'] = json.dumps(analysisMetadataDict)

	# Create canned analysis dataframe
	cannedAnalysisDataframe = mergedDataframe.merge(metadataDataframe, left_on='creeds_id', right_index=True, how='inner')[['geo_id', 'name', 'link', 'metadata']]
	cannedAnalysisDataframe = cannedAnalysisDataframe.reset_index(drop=True).rename(columns={'geo_id': 'dataset', 'name': 'tool', 'link': 'canned_analysis_url'})

	# Save
	cannedAnalysisDataframe.to_csv(outfile, sep='\t', index_label='id')

#################################################################
#################################################################
############### 2. Clustergrammer ###############################
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

##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')
