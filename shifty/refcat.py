'''
    methods to deal with Tonry's refcat2 catalog
    
    https://archive.stsci.edu/prepds/atlas-refcat2/
    http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1809.09157
'''

# -------------------------------------------------------------------------------------
# Third party imports
# -------------------------------------------------------------------------------------
import os, sys
import wget
import numpy as np
import functools
import subprocess
import pandas as pd
import io
import csv

import astropy
from astropy.wcs import WCS

# -------------------------------------------------------------------------------------
# Any local imports
# -------------------------------------------------------------------------------------
from downloader import Downloader
from pixel import Pixel

# -------------------------------------------------------------------------------------
# Various class definitions for *refcat Handling* in shifty
# -------------------------------------------------------------------------------------

class RefCat(Downloader):
    '''
        class to deal with Tonry's refcat2 catalog
        
        https://archive.stsci.edu/prepds/atlas-refcat2/
        
        methods:
        --------
        download_refcode
        download_refcat()
        compile_refcat()
        query_refcat()
        
        _read_refcat()
        
    '''
    
    def __init__(self, ) :
        # - Allow ourselves to use any Downloader methods
        Downloader.__init__(self, )
        
        # - Local directory for saving data
        self.refcat_dir = self._fetch_refcat_data_directory()
        
        # - Anticipated filepath for compiled refcat code
        self.refcat_codename = 'refcat'
        self.refcat_filepath = os.path.join(self.refcat_dir , self.refcat_codename)

        # - Set of refcat catalogues on mast
        self.downloadable_files = {
            '16':'https://archive.stsci.edu/hlsps/atlas-refcat2/orig/hlsp_atlas-refcat2_atlas_ccd_00-m-16_multi_v1_cat.tbz',
            '17':'https://archive.stsci.edu/hlsps/atlas-refcat2/orig/hlsp_atlas-refcat2_atlas_ccd_16-m-17_multi_v1_cat.tbz',
            '18':'https://archive.stsci.edu/hlsps/atlas-refcat2/orig/hlsp_atlas-refcat2_atlas_ccd_17-m-18_multi_v1_cat.tbz',
            '19':'https://archive.stsci.edu/hlsps/atlas-refcat2/orig/hlsp_atlas-refcat2_atlas_ccd_18-m-19_multi_v1_cat.tbz',
            '20':'https://archive.stsci.edu/hlsps/atlas-refcat2/orig/hlsp_atlas-refcat2_atlas_ccd_19-m-20_multi_v1_cat.tbz',
            'code' : 'https://archive.stsci.edu/prepds/atlas-refcat2/refcat.c',
            'man'  : 'https://archive.stsci.edu/prepds/atlas-refcat2/refcat.man',
        }

    # -------------------------------------------------------------------------------------
    # Public Methods
    # -------------------------------------------------------------------------------------
    def download_refcode(self, ):
        '''
           download the executable ref_cat code
        '''
        for arg in ['code', 'man']:
    
            # get the filename out of the filepath
            remote_filepath = self.downloadable_files[str(arg)]
            h,filename      = os.path.split(remote_filepath)
            
            # only bother to do the download if the file doesn't exist locally
            filepath = os.path.join(self.refcat_dir , filename )
            if not os.path.isfile(filepath):
                
                # download the file
                print(self.downloadable_files[str(arg)] , '--->>>', filepath )
                wget.download(remote_filepath, filepath)

    def download_refcat(self, *args):
        '''
            download the huge refcat catalog files
        '''
        for arg in args:
            if str(arg) in self.downloadable_files:
                print(' ------Downloading catalogue file(s) ----------')
                # get the filename out of the filepath
                remote_filepath = self.downloadable_files[str(arg)]
                h,filename      = os.path.split(remote_filepath)

                # only bother to do the download if the file doesn't exist locally
                filepath = os.path.join(self.refcat_dir , filename )
                if not os.path.isfile(filepath):
                    
                    # download the file
                    print(self.downloadable_files[str(arg)] , '--->>>', filepath )
                    wget.download(remote_filepath, filepath)
                    assert os.path.isfile(filepath)
    
                    # if a tar-file, un-tar it
                    if filename[-4:] == ".tbz":
                        self._untar_refcat( *args)

    def compile_refcat(self, ):
        '''
           try to compile the refcat code
        '''
        # change working directory to self.refcat_dir
        orig_workingdir = os.getcwd()
        os.chdir(self.refcat_dir)
        assert os.getcwd() == self.refcat_dir
        
        # instructions in refcat.man say to compile using ...
        command = 'cc -o %s -O -Wall refcat.c -lm' % self.refcat_codename
        os.system(command)
    
        # move back to original directory
        os.chdir(orig_workingdir)
        assert os.getcwd() == orig_workingdir
    
    
    def find_all_stars_on_image(self,  header, image_data , nPixels = 0 ):
        '''
            given image-data (and associated header) use refcat to find all of the stars that fall on the image
            
            inputs:
            -------
            
            returns:
            --------
            pixels / (ra,dec) / both ???
        '''

        # Identify integer RAs & Decs in the data: uses WCS:
        # Then find unique PAIRS of ra,dec integers
        sky_coords = Pixel().get_per_pixel_RADEC(header, image_data )
        int_radec = np.unique(np.array( [np.around(sky_coords.ra.degree).astype(int).flatten(),
                                         np.around(sky_coords.dec.degree).astype(int).flatten() ] ).T, axis=0)

        # search refcat around those integer ra-dec pairs
        code = self.refcat_filepath
        dir = os.path.join(self.refcat_dir, '00_m_16')
        rad = 1
        mlim = 20
        print("\
              **\
              WARNING: HARD-CODED MAG-LIMITS & SEARCH DIRECTORIES IN find_all_stars_on_image()\
              **\
              ")
        results_dict = {}
        for ra, dec in int_radec :
            results_dict.update(self._read_refcat(ra,
                                                  dec,
                                                  code,
                                                  dir,
                                                  rad=rad ,
                                                  mlim=mlim ))
        
        # get pixels corresponding to  (ra, dec)'s of catalog sources
        ra  = np.array( [v[0] for v in results_dict.values()] )
        dec = np.array( [v[1] for v in results_dict.values()] )
        pix = np.array(WCS(header).all_world2pix(ra, dec, 1))
        int_pix = np.around(pix).astype(int)
        
        # offset pixels because of offset between numpy & fits-fortran
        pix = pix-1
        int_pix = int_pix-1
        
        # remove sources which would be more than n-Pixels from the edge of the image
        ind = (int_pix[0] > -nPixels) & \
              (int_pix[0] < nPixels + header['NAXIS1']) & \
              (int_pix[1] > -nPixels) & \
              (int_pix[1] < nPixels + header['NAXIS2'])

        return  ra[ind],\
                dec[ind] ,\
                np.array( [pix[0][ind], pix[1][ind] ]),\
                np.array([int_pix[0][ind], int_pix[1][ind] ])

    # -------------------------------------------------------------------------------------
    # Private Methods
    # -------------------------------------------------------------------------------------
    def _untar_refcat(self, *args):
        '''
            un-tar a local refcat file
        '''
        for arg in args:
            if str(arg) in self.downloadable_files:
                
                # get the filename out of the filepath
                remote_filepath = self.downloadable_files[str(arg)]
                h , filename      = os.path.split(remote_filepath)
                tarred_filepath     = os.path.join(self.refcat_dir , filename )
                
                # expected un-tarred destination directory
                sub = '_ccd_'   ; start_locn = tarred_filepath.find(sub) + len(sub)
                sub = '_multi_' ; end_locn   = tarred_filepath.find(sub)
                untarred_directory = os.path.join(self.refcat_dir , tarred_filepath[start_locn:end_locn].replace('-','_') )
                
                # if the tarred file EXISTS and the UNtarred directory does NOT exist ...
                # ... then proceed to untar
                if os.path.isfile(tarred_filepath) and not os.path.isdir(untarred_directory):
                    
                    # change working directory to self.refcat_dir
                    orig_workingdir = os.getcwd()
                    os.chdir(self.refcat_dir)
                    assert os.getcwd() == self.refcat_dir

                    # instructions in refcat.man say to compile using ...
                    command = 'tar xjvf %s' % tarred_filepath
                    os.system(command)

                    # move back to original directory
                    os.chdir(orig_workingdir)
                    assert os.getcwd() == orig_workingdir
        return True

    @functools.lru_cache(maxsize=128)
    def _read_refcat(self,
                     ra,
                     dec,
                     code,
                     dir,
                     rad=3.0,
                     mlim=12.0,
                     ):
        '''
            runs refcat
            reads results from command line, each of which looks like ...
            'ra', 'dec', 'g', 'r', 'i', 'z', 'J', 'c', 'o'
            ...
            179.755775  10.937336 11.690 11.243 11.109 11.059 10.142 11.472 11.182
            ...
            
            returns:
            --------
            dictionary of results 
             - { }
        '''
        
        # run refcat and pipe results from command-line into string
        # split up string and put into dictionary
        # dictionary entries look like ...
        # key   = (str(ra),str(dec)),
        # value = np.array(['ra', 'dec', 'g', 'r', 'i', 'z', 'J', 'c', 'o'])
        result = subprocess.run([code, str(ra), str(dec),  '-rad', str(rad), '-mlim', str(mlim), '-dir', dir], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')
        return {  (row.split()[0],row.split()[1]) : np.array([float(_) for _ in row.split()  ]) for row in result if len(row)>20}
    
    

    # -------------------------------------------------------------------------------------
    # Data directories / Storage / etc
    # -------------------------------------------------------------------------------------
    def _fetch_refcat_data_directory(self):
        '''
        Create a sub-directory within self.local_dir
        Falls back to using self.local_dir
        '''

        # Define a subdirectory within the main data-directory
        refcat_filepath = os.path.join(self._fetch_data_directory(), 'refcat')

        # Attempt to create the desired sub-directory
        # Fall back to main data-directory if any error
        if not os.path.isdir(refcat_filepath):
            try:
                os.mkdir(refcat_filepath)
            except (Exception, OSError) as error :
                print(error)
                print('setting refcat download directory to be same as parent directory: %r ' % self.local_dir )
                refcat_filepath = self.local_dir
        return refcat_filepath

