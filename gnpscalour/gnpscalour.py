import webbrowser
from logging import getLogger, basicConfig, NOTSET
from logging.config import fileConfig
from pkg_resources import resource_filename
from collections import defaultdict

import numpy as np
import pandas as pd

from calour.database import Database

from . import __version_numeric__

logger = getLogger(__name__)

try:
    # get the logger config file location
    log = resource_filename(__name__, 'log.cfg')
    # set the logger output according to log.cfg
    # setting False allows other logger to print log.
    fileConfig(log, disable_existing_loggers=False)
    # set the log level to same value as calour log level
    clog = getLogger('calour')
    calour_log_level = clog.getEffectiveLevel()
    if calour_log_level != NOTSET:
        logger.setLevel(calour_log_level)
except:
    print('FAILED loading logging config file %s' % log)
    basicConfig(format='%(levelname)s:%(message)s')


class GNPS(Database):
    def __init__(self, exp=None):
        super().__init__(exp=exp, database_name='GNPS', methods=['get'])
        if 'gnpscalour' not in exp.databases:
            logger.warn('Cannot initialize GNPS database since gnps info file was not supplied. Please supply when loading')
            self.gnps_data = None
            return
        self.gnps_data = exp.databases['gnpscalour']['metabolomics_gnps_table']
        self.mzfield = 'parent mass'
        self.rtfield = 'RTMean'
        self.exp = exp
        self.gnps_data[self.mzfield] = pd.to_numeric(self.gnps_data[self.mzfield], errors='coerce')
        self.gnps_data[self.rtfield] = pd.to_numeric(self.gnps_data[self.rtfield], errors='coerce')

    def version(self):
        return __version_numeric__

    def _prepare_gnps_ids(self, direct_ids=False, mz_thresh=0.02, rt_thresh=15, use_gnps_id_from_AllFiles=True):
        '''Link each feature in the experiment to the corresponding gnps table id.

        Parameters
        ----------
        direct_ids: bool, optional
            True to link via the ids, False (default) to link via MZ/RT
        mz_thresh, rt_thresh: float, optional
            the threshold for linking to gnps if direct_ids is False
        use_gnps_id_from_AllFiles: bool, optional
            True (default) to link using the AllFiles column in GNPS, False to link using 'cluster index' column
            (if direct_ids is True).
        '''
        logger.debug('Locating GNPS ids')
        self.exp.feature_metadata['_gnps_ids'] = None
        # if we don't have the linked gnps table, all are NA
        if direct_ids:
            logger.debug('locating ids using direct_ids')
            # get the gnps ids values from the gnps file
            if use_gnps_id_from_AllFiles:
                logger.debug('using AllFiles field')
                gnps_metabolite_ids = []
                for cid in self.gnps_data['AllFiles']:
                    cid = cid.split(':')[-1]
                    cid = cid.split('###')[0]
                    cid = int(cid)
                    gnps_metabolite_ids.append(cid)
            else:
                logger.debug('using cluster index field')
                gnps_metabolite_ids = self.gnps_data['cluster index'].astype(int)
            gnps_metabolite_ids_pos = {}
            for idx, cmet in enumerate(gnps_metabolite_ids):
                gnps_metabolite_ids_pos[cmet] = idx
            gnps_ids = {}
            for cmet in self.exp.feature_metadata.index.values:
                if int(cmet) in gnps_metabolite_ids_pos:
                    gnps_ids[cmet] = [gnps_metabolite_ids_pos[int(cmet)]]
                else:
                    logger.warning('metabolite %s not found in gnps file. Are you using correct ids?' % cmet)
                    raise ValueError('metabolite ID %s not found in gnps file. Are you using correct ids?' % cmet)
        else:
            # match using MZ/RT
            logger.debug('linking using MZ (thresh %f) and RT (thresh %f)' % (mz_thresh, rt_thresh))
            gnps_ids = {}
            for cmet in self.exp.feature_metadata.index.values:
                cids = self._find_close_annotation(self.exp.feature_metadata['MZ'][cmet], self.exp.feature_metadata['RT'][cmet], mzerr=mz_thresh, rterr=rt_thresh)
                gnps_ids[cmet] = cids
        self.exp.feature_metadata['_gnps_ids'] = pd.Series(gnps_ids)

    def _prepare_gnps_names(self):
        '''Add the GNPS database name to each metabolite (if known) or 'na' if unknown.
        the metabolite name will be added in the "gnps_name" field of the feature_metadata
        the component index will be added in the "gnps_component" field of the feature_metadata
        '''
        logger.debug('Adding gnps terms as "gnps_name" and "gnps_component" columns in feature metadta')
        self.exp.add_terms_to_features('gnps', term_type='LibraryID', use_term_list=None, field_name='gnps_name')
        self.exp.add_terms_to_features('gnps', term_type='componentindex', use_term_list=None, field_name='gnps_component')
        if 'MZ' not in self.exp.feature_metadata.columns:
            logger.debug('Adding MZ field from gnps data')
            self.exp.add_terms_to_features('gnps', term_type='parent mass', use_term_list=None, field_name='MZ')
        if 'RT' not in self.exp.feature_metadata.columns:
            logger.debug('Adding RT field from gnps data')
            self.exp.add_terms_to_features('gnps', term_type='RTMean', use_term_list=None, field_name='RT')
        logger.debug('Added terms')

    def _find_close_annotation(self, mz, rt, mzerr=0.1, rterr=30):
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
        try:
            mz = float(mz)
            rt = float(rt)
            self.gnps_data.fillna(0, inplace=True)
            mzpos = np.where(np.abs(self.gnps_data[self.mzfield].values - mz) < mzerr)[0]
            rtpos = np.where(np.abs(self.gnps_data[self.rtfield] - rt) < rterr)[0]
            pos = np.intersect1d(mzpos, rtpos)
            return pos
        except:
            logger.warn('Cannot process GNPS for mz: %s, rt: %s' % (mz, rt))
            return np.empty(0)

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
            return []
        # pos = self._find_close_annotation(self._exp.feature_metadata['MZ'][feature], self._exp.feature_metadata['RT'][feature])
        pos = self._exp.feature_metadata['_gnps_ids'][feature]
        if len(pos) == 0:
            return []
        for clabel in ['parent mass', 'RTMean', 'LibraryID', 'AllOrganisms', 'componentindex']:
            shortdesc.append(({'annotationtype': 'other', 'feature': feature, 'gnps_link': self.gnps_data.iloc[pos[0]]['ProteoSAFeClusterLink']}, '%s: %s' % (clabel, self.gnps_data.iloc[pos[0]][clabel])))
        for cpos in pos:
            shortdesc.append(({'annotationtype': 'other', 'feature': feature, 'gnps_link': self.gnps_data.iloc[cpos]['ProteoSAFeClusterLink']}, str(self.gnps_data.iloc[cpos]['LibraryID'])))
        return shortdesc

    def show_annotation_info(self, annotation):
        '''Show the website for the sequence

        Parameters
        ----------
        annotation : dict
            should contain 'gnps_link'
        '''
        # open in a new tab, if possible
        new = 2

        link = self.get_annotation_website(annotation)
        if link is None:
            logger.warning('GNPS link not found. not opening.')
            return
        logger.debug('opening link: %s' % link)
        webbrowser.open(link, new=new)

    def get_annotation_website(self, annotation):
        '''Get the databasewebsite address of information about the annotation.

        Parameters
        ----------
        annotation : dict
            keys/values are database specific.
            E.g. See dbBact REST API /annotations/get_annotation for keys / values


        Returns
        -------
        str or None
            The webaddress of the html page with details about the annotation,
            or None if not available
        '''
        if 'gnps_link' not in annotation:
            logger.warning('GNPS info does not contain gnps_link. It only contains: %s' % list(annotation.keys()))
            return None
        link = annotation['gnps_link'] + '&show=True'
        return link

    def get_feature_terms(self, features, exp=None, term_type='LibraryID', ignore_exp=None):
        '''Get list of gnps terms per feature.

        Parameters
        ----------
        features : list of str
            the features to get the terms for
        exp : calour.Experiment (optional)
            not None to store results inthe exp (to save time for multiple queries)
        term_type : str, optional
            the GNPS cluster file column names. common options include:
                'LibraryID': (default) the metabolite name
                'componentindex': the cluster to which this metabolite belongs
                'AllOrganisms': the source library?
        ignore_exp : not used

        Returns
        -------
        feature_terms: dict of dict of float
            key is the feature, value is dict of key: gnps id, value is 1 (to be compatible with calour.database)
        '''
        if term_type not in self.gnps_data.columns:
            raise ValueError('term_type %s not a column in the gnps clusterinfo file.' % term_type)
        feature_terms = defaultdict(list)
        for cfeature in features:
            # pos = self._find_close_annotation(self._exp.feature_metadata['MZ'][cfeature], self._exp.feature_metadata['RT'][cfeature])
            pos = self._exp.feature_metadata['_gnps_ids'][cfeature]
            cterms = {}
            foundna = False
            for cpos in pos:
                cterm = self.gnps_data.iloc[cpos][term_type]
                # pandas treats N/A as nan
                # if cterm == 'N/A' or not isinstance(cterm, str):
                if cterm == 'N/A':
                    foundna = True
                    continue
                cterms[cterm] = 1
            feature_terms[cfeature] = cterms
            if len(cterms) == 0:
                if foundna:
                    feature_terms[cfeature] = ['na']
        return feature_terms
