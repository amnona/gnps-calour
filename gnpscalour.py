import webbrowser
from logging import getLogger

import numpy as np

from calour.database import Database

logger = getLogger(__name__)


class GNPS(Database):
    def __init__(self, exp=None):
        super().__init__(exp=exp, database_name='GNPS', methods=['get'])
        if '_calour_metabolomics_gnps_table' not in exp.exp_metadata:
            raise ValueError('Cannot initialize GNPS database since gnps info file was not supplied. Please supply when loading ')
        self.gnps_data = exp.exp_metadata['_calour_metabolomics_gnps_table']
        self.mzfield = 'parent mass'
        self.rtfield = 'RTMean'

    def _find_close_annotation(self, mz, rt, mzerr=20, rterr=100):
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

        address = annotation['gnps_link']
        webbrowser.open(address, new=new)
