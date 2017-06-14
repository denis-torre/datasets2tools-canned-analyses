import gzip, urllib2, os
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

            # Get platform
            self.platform_accession = os.path.basename(matrix_url).split('_')[1]

#######################################################
########## 2. Get Expression ##########################
#######################################################

    def get_expression_data(self, n_genes=500):

        # Check platform
        if self.platform_accession in ['GPL11154', 'GPL13112', 'GPL16791', 'GPL17021']:

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

    def get_sample_metadata(self):

        # Get metadata dataframe
        metadata_dataframe = pd.DataFrame([x.split('\t') for x in self.data.split('\n')[1:] if any(y in x for y in ['!Sample', '!^SAMPLE'])]).set_index(0)

        # Get title conversion
        sample_title_dict = {sample_accession: '{sample_title} ({sample_accession})'.format(**locals()) for sample_accession, sample_title in metadata_dataframe.loc[['!^SAMPLE', '!Sample_title']].T.as_matrix() if sample_accession}

        # Get metadata dict
        metadata_dict = [{term_string.split(': ')[0]: term_string.split(': ')[1] for term_string in term_list} for term_list in np.array_split(metadata_dataframe.loc['!Sample_characteristics_ch1'].dropna().tolist(), len(sample_title_dict.keys()))]

        # Create dict
        sample_metadata_dataframe = pd.DataFrame({sample_title_dict[sample_accession]: metadata_dict for sample_accession, metadata_dict in zip(metadata_dataframe.loc['!^SAMPLE'], metadata_dict)})

        # Rename
        self.expression_dataframe = self.expression_dataframe.rename(columns=sample_title_dict)

        # Add cats
        self.sample_cats = [{'title': index, 'cats': {value: rowData[rowData==value].index.tolist() for value in set(rowData.values)}} for index, rowData in sample_metadata_dataframe.iterrows()]

#######################################################
########## 4. Get JSON ################################
#######################################################

    def get_clustergrammer_json(self, outfile):

        # Create network
        net = Network()

        # Load file
        net.load_df(self.expression_dataframe)

        # Add categories
        try:
            net.add_cats('col', self.sample_cats)
        except:
            pass

        try:
            # calculate clustering using default parameters
            net.cluster()

            # save visualization JSON to file for use by front end
            net.write_json_to_file('viz', outfile)
        except:
            os.system('touch {outfile}'.format(**locals()))


#######################################################
########## .  ##################################
#######################################################
