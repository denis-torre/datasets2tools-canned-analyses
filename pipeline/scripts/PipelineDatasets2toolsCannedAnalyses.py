#!/usr/bin/env python
# -*- coding: utf-8 -*- 
#################################################################
#################################################################
###############  ################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
import sys, os
import pandas as pd
from clustergrammer_widget import *
import numpy as np
import requests

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
import Support as S 

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####

#######################################################
#######################################################
########## S1. CREEDS
#######################################################
#######################################################

#############################################
########## 1. Get Analysis Title
#############################################

def getAnalysisTitle(toolName, metadataDict):
    
    # Tool
    if toolName == 'enrichr':
        title = 'Enrichment analysis of genes ' + metadataDict['geneset']
        if 'disease_name' in metadataDict.keys():
            title += ' in ' + metadataDict['disease_name']
        elif 'drug_name' in metadataDict.keys():
            title += ' following ' + metadataDict['drug_name'] + ' treatment'
        elif 'hs_gene_symbol' in metadataDict.keys():
            perturbationLabel = metadataDict['pert_type'] if 'pert_type' in metadataDict.keys() else 'perturbation'
            title += ' in  ' + metadataDict['hs_gene_symbol'] + ' ' + perturbationLabel
    elif toolName == 'paea':
        title = 'Enrichment analysis of genes dysregulated '
        if 'disease_name' in metadataDict.keys():
            title += ' in ' + metadataDict['disease_name']
        elif 'drug_name' in metadataDict.keys():
            title += ' following ' + metadataDict['drug_name'] + ' treatment'
        elif 'hs_gene_symbol' in metadataDict.keys():
            perturbationLabel = metadataDict['pert_type'] if 'pert_type' in metadataDict.keys() else 'perturbation'
            title += ' in  ' + metadataDict['hs_gene_symbol'] + ' ' + perturbationLabel
    elif toolName == 'l1000cds2':
        title = 'Small molecules which ' + metadataDict['direction'] + ' '
        if 'disease_name' in metadataDict.keys():
            title += metadataDict['disease_name']
        elif 'drug_name' in metadataDict.keys():
            title += metadataDict['drug_name'] + ' treatment'
        elif 'hs_gene_symbol' in metadataDict.keys():
            perturbationLabel = metadataDict['pert_type'] if 'pert_type' in metadataDict.keys() else 'perturbation'
            title += metadataDict['hs_gene_symbol'] + ' ' + perturbationLabel
    elif toolName == 'genemania':
        title = 'Interaction network and enrichment analysis of genes ' + metadataDict['geneset']
        if 'disease_name' in metadataDict.keys():
            title += ' in ' + metadataDict['disease_name']
        elif 'drug_name' in metadataDict.keys():
            title += ' following ' + metadataDict['drug_name'] + ' treatment'
        elif 'hs_gene_symbol' in metadataDict.keys():
            perturbationLabel = metadataDict['pert_type'] if 'pert_type' in metadataDict.keys() else 'perturbation'
            title += ' in ' + metadataDict['hs_gene_symbol'] + ' ' + perturbationLabel

    return title

#############################################
########## 2. Get Analysis Description
#############################################

def getAnalysisDescription(toolName, metadataDict):

    # Cell type
    celltypeLabel = ' '.join([' in the ', metadataDict['cell_type'], ' cell type.']) if 'cell_type' in metadataDict.keys() else '.'

    # Tool
    if toolName == 'enrichr':
        title = 'An enrichment analysis was performed on the top 500 most ' + metadataDict['geneset'] + ' genes identified by applying the Characteristic Direction method comparing gene expression between '
        if 'disease_name' in metadataDict.keys():
            title += 'cells affected by ' + metadataDict['disease_name'] + ' and healthy control cells' + celltypeLabel
        elif 'drug_name' in metadataDict.keys():
            title += metadataDict['drug_name'] + ' treated cells and untreated control cells' + celltypeLabel
        elif 'hs_gene_symbol' in metadataDict.keys():
            perturbationLabel = metadataDict['pert_type'] if 'pert_type' in metadataDict.keys() else 'perturbation'
            title += metadataDict['hs_gene_symbol'] + ' ' + perturbationLabel + ' cells and control cells' + celltypeLabel
    elif toolName == 'paea':
        title = 'An enrichment analysis was performed on the top most dyresgulated genes determined by applying the principal angle method to compare gene expression between '
        if 'disease_name' in metadataDict.keys():
            title += 'cells affected by ' + metadataDict['disease_name'] + ' and healthy control cells' + celltypeLabel
        elif 'drug_name' in metadataDict.keys():
            title += metadataDict['drug_name'] + ' treated cells and untreated control cells' + celltypeLabel
        elif 'hs_gene_symbol' in metadataDict.keys():
            perturbationLabel = metadataDict['pert_type'] if 'pert_type' in metadataDict.keys() else 'perturbation'
            title += metadataDict['hs_gene_symbol'] + ' ' + perturbationLabel + ' cells and control cells' + celltypeLabel
    elif toolName == 'l1000cds2':
        title = 'The L1000 database was queried in order to identify small molecule perturbations which ' + metadataDict['direction'] + ' the '
        if 'disease_name' in metadataDict.keys():
            title += metadataDict['disease_name'] + ' gene expression signature' + celltypeLabel
        elif 'drug_name' in metadataDict.keys():
            title += ' gene expression signature induced by ' + metadataDict['drug_name'] + ' treatment' + celltypeLabel
        elif 'hs_gene_symbol' in metadataDict.keys():
            perturbationLabel = metadataDict['pert_type'] if 'pert_type' in metadataDict.keys() else 'perturbation'
            title += ' gene expression signature induced by ' + metadataDict['hs_gene_symbol'] + ' ' + perturbationLabel + celltypeLabel
    elif toolName == 'genemania':
        title = 'The analysis explores the gene interaction network and pathway enrichment of the top 50 most ' + metadataDict['geneset'] + ' genes identified by applying the Characteristic Direction method comparing gene expression between '
        if 'disease_name' in metadataDict.keys():
            title += 'cells affected by ' + metadataDict['disease_name'] + ' and healthy control cells' + celltypeLabel
        elif 'drug_name' in metadataDict.keys():
            title += metadataDict['drug_name'] + ' treated cells and untreated control cells' + celltypeLabel
        elif 'hs_gene_symbol' in metadataDict.keys():
            perturbationLabel = metadataDict['pert_type'] if 'pert_type' in metadataDict.keys() else 'perturbation'
            title += metadataDict['hs_gene_symbol'] + ' ' + perturbationLabel + ' cells and control cells' + celltypeLabel

    return title

#############################################
########## 1. Get CREEDS Description
#############################################

def generateDescriptions(cannedAnalysisDataframe):
    metadata_list = []
    for tool_name, metadata in cannedAnalysisDataframe[['tool_name', 'metadata']].as_matrix():
        if tool_name == 'enrichr' or tool_name == 'paea':
            if 'disease_name' in metadata.keys():
                enrichment_signature_label = ' in ' + metadata['disease_name']
            elif 'drug_name' in metadata.keys():
                enrichment_signature_label = ' following ' + metadata['drug_name'] + ' treatment'
            elif 'hs_gene_symbol' in metadata.keys() and 'pert_type' in metadata.keys():
                enrichment_signature_label = ' following ' + ' '.join([metadata['hs_gene_symbol'], metadata['pert_type']]) 
            else:
                enrichment_signature_label = ''
        else:
            if 'disease_name' in metadata.keys():
                l1000cds2_signature_label = metadata['disease_name']
            elif 'drug_name' in metadata.keys():
                l1000cds2_signature_label = metadata['drug_name'] + ' treatment'
            elif 'hs_gene_symbol' in metadata.keys() and 'pert_type' in metadata.keys():
                l1000cds2_signature_label = ' '.join([metadata['hs_gene_symbol'], metadata['pert_type']]) 
            else:
                l1000cds2_signature_label = ''
        if 'cell_type' in metadata.keys():
            cell_type_label = ' in ' + metadata['cell_type']
        if tool_name == 'enrichr':
            metadata['description'] = 'Enrichment analysis of genes ' + metadata['geneset'] + enrichment_signature_label + cell_type_label
        elif tool_name == 'l1000cds2':
            metadata['description'] = 'Query for small molecules which ' + ' '.join([metadata['direction'], l1000cds2_signature_label]) + ' gene expression signature' + cell_type_label
        elif tool_name == 'paea':
            metadata['description'] = 'Enrichment analysis of differential gene expression signature' + enrichment_signature_label + cell_type_label
        metadata_list.append(json.dumps(metadata))
    return metadata_list

#######################################################
#######################################################
########## S2. GEO Matrix
#######################################################
#######################################################

#############################################
########## 1. Get Matrix Link
#############################################

def getMatrixLink(geo_accession, gpl, ftp_link):
    geo_type = geo_accession[:3]
    gpl_list = ['GPL'+x for x in gpl.split(';')]
    if geo_type == 'GDS':
        result_list = [[geo_accession, gpl_list[0], os.path.join(ftp_link, 'soft', geo_accession+'_full.soft.gz')]]
    elif geo_type == 'GSE':
        if len(gpl_list) > 1:
            result_list = [[geo_accession, gpl_str, os.path.join(ftp_link, 'matrix', geo_accession+'-'+gpl_str+'_series_matrix.txt.gz')] for gpl_str in gpl_list]
        else:
            result_list = [[geo_accession, gpl_list[0], os.path.join(ftp_link, 'matrix', geo_accession+'_series_matrix.txt.gz')]]
    else:
        raise ValueError('GEO type incorrect (must be GDS or GSE).')
    return result_list

#######################################################
#######################################################
########## S3. GEO Clustergram - Soft
#######################################################
#######################################################

#############################################
########## 1. Soft File
#############################################

def prepareSoftFile(geoData, nGenes=1000):
    
    # Define index
    index = 'IDENTIFIER'
    
    # Get matrix
    splitData = [x.split('\t') for x in geoData if x[0] not in ['^', '!', '#']]
    expressionData = pd.DataFrame(splitData[1:], columns=splitData[0])

    # Replace and drop NA
    expressionData = expressionData.replace('null', np.nan).replace('NULL', np.nan).replace('', np.nan).dropna()

    # Fix columns
    columnsToKeep = [x for x in expressionData.columns if x[:3] == 'GSM' or x == index]
    expressionData = expressionData[columnsToKeep]
    
    # Group by and aggregate by mean
    expressionData = expressionData.set_index(index).astype('float').reset_index().groupby(index).mean()

    # Get sorted genes by variance
    topGenes = expressionData.apply(np.var, 1).sort_values(ascending=False).index[:nGenes]
    
    # Get subset
    expressionData = expressionData.loc[topGenes]
    
    # Fix index
    expressionData.index = ['Gene Symbol: '+x.split('//')[0].replace(' ','') for x in expressionData.index]
       
    # Get sample titles
    sampleTitleDict = {y.split(':')[0]: y.replace(':', '׃') for y in [x.split('Value for ')[1] for x in geoData if x[:4] == '#GSM']}

    # Define annotation dict
    annotationList = []

    # Loop through lines
    for geoLine in geoData:

        # Get subset description
        if geoLine[:19] == '!subset_description':

            # Get description string
            subsetDescription = geoLine.split('= ')[1]

        # Get subset samples
        if geoLine[:17] == '!subset_sample_id':

            # Get sample string
            subsetSamples = geoLine.split('= ')[1].split(',')

            # Add 
            annotationList += [[subsetDescription, sample, 'true'] for sample in subsetSamples]

    # Convert to dataframe
    annotationDataframe = pd.DataFrame(annotationList).drop_duplicates().pivot(index=0, columns=1, values=2).fillna('false')

    # Convert to dict
    columnAnnotationDict = annotationDataframe.to_dict()

    # Convert to tuples
    columnTupleDict = {x:tuple(['Sample: ' + sampleTitleDict[x]]+[': '.join([key, str(item)]) for key, item in columnAnnotationDict[x].iteritems()]) for x in columnAnnotationDict.keys()}
    
    # Fix colnames
    expressionData.columns = [columnTupleDict[x] for x in expressionData.columns]
    
    # Return
    return expressionData

#############################################
########## 2. Series Matrix
#############################################

def prepareSeriesMatrix(geoData, platformFile, nGenes=1000):
    
    # Split and replace
    geoDataProcessed = [x.replace('"', '').split('\t') for x in geoData]

    # Get expression data
    expressionData = [x for x in geoDataProcessed if len(x) > 1 and x[0][0] != '!']

    # Make expression dataframe
    expressionDataframe = pd.DataFrame(expressionData[1:], columns=expressionData[0])

    # Read platform dataframe
    platformDataframe = pd.read_table(platformFile, comment='#')#, usecols=['ID', symbolColumn]).replace('---',np.nan).dropna()

    # Get symbol columns
    platform = os.path.basename(platformFile).split('-')[0]
    if platform in ['GPL10558', 'GPL6885', 'GPL6887', 'GPL6947']:
        symbolColumn = 'Symbol'
    elif platform in ['GPL96', 'GPL570', 'GPL1261', 'GPL339', 'GPL570', 'GPL571', 'GPL81', 'GPL8321']:
        symbolColumn = 'Gene Symbol'
    elif platform in ['GPL4133', 'GPL4134', 'GPL7202', 'GPL6480']:
        symbolColumn = 'GENE_SYMBOL'
    elif platform in ['GPL6244', 'GPL6246']:
        platformDataframe['Symbol'] = [x.split(' // ')[1] if type(x) == str and ' // ' in x else None for x in platformDataframe['gene_assignment']]
        symbolColumn = 'Symbol'

    # Filter
    platformDataframe = platformDataframe[['ID', symbolColumn]].dropna()

    # Fix
    platformDataframe['ID'] = [str(x) for x in platformDataframe['ID']]
    expressionDataframe['ID_REF'] = [str(x) for x in expressionDataframe['ID_REF']]
    expressionDataframe = expressionDataframe.replace('', np.nan).replace('null', np.nan).replace('NULL', np.nan).dropna()

    # Merge
    expressionDataframe = platformDataframe.merge(expressionDataframe, left_on='ID', right_on='ID_REF', how='inner').drop(['ID', 'ID_REF'], axis=1)

    # Collapse
    expressionDataframe = expressionDataframe.set_index(symbolColumn).astype('float').reset_index().groupby(symbolColumn).mean()

    # Get sorted genes by variance
    topGenes = expressionDataframe.apply(np.var, 1).sort_values(ascending=False).index[:nGenes]

    # Get subset
    expressionDataframe = expressionDataframe.loc[topGenes]

    # Fix index
    expressionDataframe.index = ['Gene Symbol: '+x.split('//')[0].replace(' ','') for x in expressionDataframe.index]

    # Get selected columns
    selectedColumns = ['!Sample_geo_accession', '!Sample_source_name_ch1', '!Sample_organism_ch1', '!Sample_characteristics_ch1']

    # Convert to dataframe
    annotationDataframe = pd.DataFrame([x for x in geoDataProcessed if x[0] in selectedColumns]).set_index(0)

    # Fix columns
    annotationDataframe.columns = [accession+'׃ '+name for accession, name in annotationDataframe.loc[['!Sample_geo_accession', '!Sample_source_name_ch1'],:].T.as_matrix()]
    annotationDataframe.drop(['!Sample_geo_accession', '!Sample_source_name_ch1'], inplace=True)

    # Fix organism
    if '!Sample_organism_ch1' in annotationDataframe.index:
        annotationDataframe.loc['!Sample_organism_ch1'] = ['organism: '+x for x in annotationDataframe.loc['!Sample_organism_ch1']]

    # Drop non variable rows and rows without a single :
    annotationCounts = annotationDataframe.apply(lambda x: len(set(x)), 1)
    annotationDataframeFiltered = annotationDataframe.iloc[[i for i, e in enumerate(annotationCounts) if e > 1]]
    annotationDataframeFiltered.drop([index for index, rowData in annotationDataframeFiltered.iterrows() if not any([':' in x for x in rowData])], inplace=True)

    if len(annotationDataframeFiltered.index) > 0:
    # Melt
        annotationDataframeMelt = pd.melt(annotationDataframeFiltered, var_name='sample_title', value_name='sample_characteristics').replace('', np.nan).dropna()

        # Get values
        annotationDataframeMelt['variable'] = [x.split(':', 1)[0] for x in annotationDataframeMelt['sample_characteristics']]
        annotationDataframeMelt['value'] = [x.split(':', 1)[1] for x in annotationDataframeMelt['sample_characteristics']]

        # Cast
        annotationDataframeCast = annotationDataframeMelt.drop_duplicates(['sample_title', 'variable']).pivot(index='variable', columns='sample_title', values='value')

        # Fill empty cells
        annotationDataframeCast = annotationDataframeFiltered.T.merge(annotationDataframeCast.T, left_index=True, right_index=True, how='left').T.drop(annotationDataframeFiltered.index).fillna('false')

        # Convert to dict
        columnAnnotationDict = annotationDataframeCast.fillna('false').to_dict()

        # Convert to tuples
        columnTupleDict = {x.split('׃')[0]:tuple(['Sample: ' + x]+[': '.join([key, str(item)]) for key, item in columnAnnotationDict[x].iteritems()]) for x in columnAnnotationDict.keys()}

        # Add columns which aren't specified
        for column in annotationDataframeFiltered.columns:
            sampleId = column.split('׃')[0]
            if sampleId not in columnTupleDict.keys():
                columnTupleDict[sampleId] = column

    else:
        columnTupleDict = {x.split('׃')[0]:'Sample: ' + x for x in annotationDataframeFiltered.columns}

    # Fix colnames
    expressionDataframe.columns = [columnTupleDict[x] for x in expressionDataframe.columns]

    # Return
    return expressionDataframe

#############################################
########## 3. Get Clustergrammer Link
#############################################

def getClustergramLink(matrixFile):

    # Read data
    with open(matrixFile, 'r') as openfile:
        geoData = [x.strip() for x in openfile.readlines()]

    # Get expression data
    expressionData = getMatrix(geoData)

    # Get annotation data
    columnAnnotationDict = getColumnAnnotations(geoData)

    # Fix columns
    expressionData.columns = [columnAnnotationDict[x] for x in expressionData.columns]

    # Get filename
    tempfile = os.path.basename(matrixFile).split('_')[0]+'.txt'

    # Save file
    expressionData.to_csv(tempfile, sep='\t', index_label='Gene Symbol')

    # Get link
    upload_url = 'http://amp.pharm.mssm.edu/clustergrammer/matrix_upload/'

    # Post request
    r = requests.post(upload_url, files={'file': open(tempfile, 'rb')})

    # Get link string
    link = os.path.dirname(r.text)

    # Remove file
    os.unlink(tempfile)

    # Return
    return link

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################

