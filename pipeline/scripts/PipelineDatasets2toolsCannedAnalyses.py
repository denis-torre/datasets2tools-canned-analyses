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
########## S1. GEO
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
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################

