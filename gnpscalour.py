import webbrowser
from logging import getLogger
from collections import defaultdict

import numpy as np

from calour.database import Database

logger = getLogger(__name__)


class GNPS(Database):
    def __init__(self, exp=None):
        super().__init__(exp=exp, database_name='GNPS', methods=['get'])
        if '_calour_metabolomics_gnps_table' not in exp.exp_metadata:
            logger.warn('Cannot initialize GNPS database since gnps info file was not supplied. Please supply when loading')
            self.gnps_data = None
            return
        self.gnps_data = exp.exp_metadata['_calour_metabolomics_gnps_table']
        self.mzfield = 'parent mass'
        self.rtfield = 'RTMean'

    def _find_close_annotation(self, mz, rt, mzerr=0.2, rterr=10):
        '''Find gnps annotations with mz,rt close enough to the requested mz,rt

        Parameters
        ----------
        mz : float
            the m/z we're looking for
        rt : float
            the retention time we're looking for
        mzerr : float (optional)
            the error window for mz
        rterr : float (optional)
            the error window for rt

        Returns
        -------
        list of close annotation indices
        '''
        mzpos = np.where(np.abs(self.gnps_data[self.mzfield].values - mz) < mzerr)[0]
        rtpos = np.where(np.abs(self.gnps_data[self.rtfield] - rt) < rterr)[0]
        pos = np.intersect1d(mzpos, rtpos)
        return pos

    def get_seq_annotation_strings(self, feature):
        '''Get nice string summaries of annotations for a given sequence

        Parameters
        ----------
        sequence : str
            the DNA sequence to query the annotation strings about

        Returns
        -------
        shortdesc : list of (dict,str) (annotationdetails,annotationsummary)
            a list of:
                annotationdetails : dict
                    'seqid' : str, the sequence annotated
                    'annotationtype : str
                    ...
                annotationsummary : str
                    a short summary of the annotation
        '''
        shortdesc = []
        if self.gnps_data is None:
            shortdesc = []
            return
        pos = self._find_close_annotation(self._exp.feature_metadata['MZ'][feature], self._exp.feature_metadata['RT'][feature])
        for cpos in pos:
            shortdesc.append(({'annotationtype': 'other', 'feature': feature, 'gnps_link': self.gnps_data.iloc[cpos]['ProteoSAFeClusterLink']}, str(self.gnps_data.iloc[cpos]['LibraryID'])))
        return shortdesc

    def show_annotation_info(self, annotation):
        '''Show the website for the sequence

        Parameters
        ----------
        annotation : dict
            should contain 'sequence'
        '''
        # open in a new tab, if possible
        new = 2

        address = annotation['gnps_link']+'&show=True'
        webbrowser.open(address, new=new)

    def get_feature_terms(self, features, exp=None, term_type=None):
        '''Get list of gnps terms per feature.

        Parameters
        ----------
        features : list of str
            the features to get the terms for
        exp : calour.Experiment (optional)
            not None to store results inthe exp (to save time for multiple queries)
        term_type : not used

        Returns
        -------
        feature_terms : dict of list of str/int
            key is the feature, list contains all terms associated with the feature
        '''
        feature_terms = defaultdict(list)
        for cfeature in features:
            pos = self._find_close_annotation(self._exp.feature_metadata['MZ'][cfeature], self._exp.feature_metadata['RT'][cfeature])
            cterms = []
            foundna = False
            for cpos in pos:
                cterm = self.gnps_data.iloc[cpos]['LibraryID']
                # pandas treats N/A as nan
                if cterm == 'N/A' or not isinstance(cterm, str):
                    foundna = True
                    continue
                cterms.append(cterm)
            feature_terms[cfeature] = cterms
            if len(cterms) == 0:
                if foundna:
                    feature_terms[cfeature] = ['na']
        return feature_terms
