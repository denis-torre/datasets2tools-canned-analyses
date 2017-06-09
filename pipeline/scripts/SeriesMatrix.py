import gzip, urllib2
import pandas as pd
import numpy as np
from StringIO import StringIO
from clustergrammer import Network

class SeriesMatrix:

#######################################################
########## 1. Initialize ##############################
#######################################################

    def __init__(self, matrix_url):

            # Open URL
            request = urllib2.Request(matrix_url)

            # Add header
            request.add_header('Accept-encoding', 'gzip')

            # Read response
            response = urllib2.urlopen(request)

            # Convert response
            buf = StringIO(response.read())

            # Open gzip file
            f = gzip.GzipFile(fileobj=buf)

            # Read data
            self.data = f.read()

#######################################################
########## 2. Get Expression ##########################
#######################################################

    def get_expression_data(self, n_genes=500):

        # Get expression dataframe
        expression_dataframe = pd.DataFrame([x.split('\t') for x in self.data.split('\n')[1:] if '!' not in x])

        # Fix axis names
        expression_dataframe = expression_dataframe.rename(columns=expression_dataframe.iloc[0]).drop(0).set_index('ID_REF').fillna(0).astype('int')

        # Get variable genes
        top_variance_genes = expression_dataframe.apply(np.var, 1).sort_values(ascending=False).index.tolist()[:n_genes]

        # Filter dataframe
        expression_dataframe = expression_dataframe.loc[top_variance_genes]

        # Z-score
        self.expression_dataframe = ((expression_dataframe.T - expression_dataframe.T.mean())/expression_dataframe.T.std()).T

#######################################################
########## 3. Get Metadata ############################
#######################################################

    def get_sample_metadata(self, terms=['!Sample_title', '!Sample_source_name_ch1', '!Sample_characteristics_ch1']):

        # Get metadata dataframe
        metadata_dataframe = pd.DataFrame([x.split('\t') for x in self.data.split('\n')[1:] if any(y in x for y in ['!Sample', '!^SAMPLE'])]).set_index(0)

        # Fix column names
        metadata_dataframe = metadata_dataframe.rename(columns=metadata_dataframe.loc['!^SAMPLE']).drop('!^SAMPLE').drop(None, axis=1)

        # Get rows of interest
        metadata_rows = set(metadata_dataframe.index).intersection(set(terms))

        # Get subset
        metadata_dataframe = metadata_dataframe.loc[metadata_rows, self.expression_dataframe.columns.tolist()]

        # Cats list
        cats_list = []

        # Loop through dataframe
        for catName, rowData in metadata_dataframe.iterrows():
            
            # Define dict
            sample_cats_dict = {}
          
            # Reindex
            metadata_dataframe_reindexed = metadata_dataframe.loc[catName].to_frame().reset_index().set_index(catName)
             
            # Get unique category values
            catValues = set(rowData.values)
             
            # Loop through values
            for catValue in catValues:
                
                # Get sample accessions
                sample_accessions = metadata_dataframe_reindexed.loc[catValue, 'index']
                
                # Convert to list
                sample_accessions = [sample_accessions] if type(sample_accessions) == str else list(sample_accessions)
                
                # Add
                sample_cats_dict[catValue] = sample_accessions
                
            # Append
            cats_list.append({'title': catName, 'cats': sample_cats_dict})

        # Add to self
        self.sample_cats = cats_list

#######################################################
########## 4. Get Metadata ############################
#######################################################

    def get_clustergrammer_json(self, outfile):

        # Create network
        net = Network()

        # Load file
        net.load_df(self.expression_dataframe)

        # Add categories
        net.add_cats('col', self.sample_cats)

        # calculate clustering using default parameters
        net.cluster()

        # save visualization JSON to file for use by front end
        net.write_json_to_file('viz', outfile)


#######################################################
########## .  ##################################
#######################################################
