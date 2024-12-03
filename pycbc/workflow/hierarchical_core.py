# Copyright (C) 2013, 2017  Ian Harry, Alex Nitz, Duncan Brown
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
This module provides the worker functions and classes that are used when
creating a workflow. For details about the workflow module see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope.html
"""
import os, stat, subprocess, logging, math, string, urllib, pickle, copy
import configparser as ConfigParser
from urllib.request import pathname2url
from urllib.parse import urljoin
import numpy, random
from itertools import combinations, groupby, permutations
from operator import attrgetter
import lal
import lal.utils
import Pegasus.api  # Try and move this into pegasus_workflow
from glue import lal as gluelal
from ligo import segments
from ligo.lw import lsctables, ligolw
from ligo.lw import utils as ligolw_utils
from ligo.lw.utils import segments as ligolw_segments
from pycbc import makedir
from pycbc.io.ligolw import LIGOLWContentHandler, create_process_table
from pycbc.workflow import pegasus_workflow
from pycbc.workflow.configuration import WorkflowConfigParser, resolve_url
from pycbc.workflow.pegasus_sites import make_catalog


def make_analysis_dir(path):
    """
    Make the analysis directory path, any parent directories that don't already
    exist, and the 'logs' subdirectory of path.
    """
    if path is not None:
        makedir(os.path.join(path, 'logs'))
        

file_input_from_config_dict = {}

# Kanchan: Executable is modified for h-search (pycbc: 2.0.1? commit 0.0a8787)
class Executable(pegasus_workflow.Executable):
    # These are the file retention levels
    DO_NOT_KEEP = 0
    INTERMEDIATE_PRODUCT = 1
    ALL_TRIGGERS = 2
    MERGED_TRIGGERS = 3
    FINAL_RESULT = 4

    # Set this parameter to indicate that this option is used to specify a
    # file and is *not* handled explicitly in the create_node or __init__
    # methods of the sub-class. Usually that is to say that this option is a
    # file and is normally specified in an file, e.g. a PSD file. As files
    # need to be identified as such to pegasus, this attempts to catch this
    # case.
    # These are standard file input arguments used in PyCBC, so we declare
    # these as files if given to any PyCBC job.
    file_input_options = ['--gating-file', '--frame-files', '--injection-file',
                          '--statistic-files', '--bank-file', '--config-files',
                          '--psd-file', '--asd-file',
                          '--fake-strain-from-file',
                          '--sgburst-injection-file']


    # Set this parameter to indicate that this option should take different
    # values based on the time. E.g. something like
    # --option1 value1[0:1000],value2[1000:2000]
    # would be replaced with --option1 value1 if the time is within 0,1000 and
    # value2 if in 1000,2000. A failure will be replaced if the job time is
    # not fully contained in one of these windows, or if fully contained in
    # multiple of these windows. This is resolved when creating the Job from
    # the Executable
    time_dependent_options = []

    # This is the default value. It will give a warning if a class is
    # used where the retention level is not set. The file will still be stored
    KEEP_BUT_RAISE_WARNING = 5
    _warned_classes_list = ['Executable']

    # Sub classes, or instances, should override this. If not overriden the
    # file will be retained, but a warning given
    current_retention_level = KEEP_BUT_RAISE_WARNING
    def __init__(self, cp, name, ifos=None, out_dir=None, tags=None,
                 reuse_executable=True, set_submit_subdir=True):
        """
        Initialize the Executable class.

        Parameters
        -----------
        cp : ConfigParser object
            The ConfigParser object holding the workflow configuration settings
        exec_name : string
            Executable name
        universe : string, optional
            Condor universe to run the job in
        ifos : string or list, optional
            The ifo(s) that the Job is valid for. If the job is
            independently valid for multiple ifos it can be provided as a list.
            Ie. ['H1',L1','V1'], if the job is only valid for the combination
            of ifos (for e.g. ligolw_thinca) then this can be supplied
            as, for e.g. "H1L1V1".
        out_dir: path, optional
            The folder to store output files of this job.
        tags : list of strings
            A list of strings that is used to identify this job.
        """
        if isinstance(ifos, str):
            self.ifo_list = [ifos]
        else:
            self.ifo_list = ifos
        if self.ifo_list is not None:
            self.ifo_list = sorted(self.ifo_list)
            self.ifo_string = ''.join(self.ifo_list)
        else:
            self.ifo_string = None
        self.cp = cp
        self.name = name
        self.container_cls = None
        self.container_type = None

        try:
            self.installed = cp.getboolean('pegasus_profile-%s' % name,
                                           'pycbc|installed')
        except:
            self.installed = False

        self.update_current_tags(tags)

        self.update_output_directory(out_dir=out_dir)

        # Determine the level at which output files should be kept
        if cp.has_option_tags('pegasus_profile-%s' % name,
                              'pycbc|retention_level', tags):
            # Get the retention_level from the config file
            # This method allows us to use the retention levels
            # defined above
            cfg_ret_level = cp.get_opt_tags('pegasus_profile-%s' % name,
                    'pycbc|retention_level', tags)
            self.current_retention_level = getattr(self, cfg_ret_level)

        self.update_current_retention_level(self.current_retention_level)

        # Should I reuse this executable?
        if reuse_executable:
            self.pegasus_name = self.name
        else:
            self.pegasus_name = self.tagged_name

        # Check that the executable actually exists locally or
        # looks like a URL, in which case trust Pegasus to be
        # able to fetch it.
        exe_path = cp.get('executables', name)
        self.needs_fetching = False

        exe_url = urllib.parse.urlparse(exe_path)

        # See if the user specified a list of sites for the executable
        # Ordering is:
        #  1) Check if a specific site for this Executable is set.
        #  2) Check is primary_site is set globally.
        #  3) Use condorpool_symlink as a fallback.
        self.exe_pfns = {}
        if cp.has_option_tags('pegasus_profile-%s' % name, 'pycbc|site', tags):
            exe_site = cp.get_opt_tags('pegasus_profile-%s' % name,
                                       'pycbc|site', tags)
        elif cp.has_option('pegasus_profile', 'pycbc|primary_site'):
            exe_site = cp.get('pegasus_profile', 'pycbc|primary_site')
        else:
            exe_site = 'condorpool_symlink'

        exe_site = exe_site.strip()

        if exe_url.scheme in ['', 'file']:
            # NOTE: There could be a case where the exe is available at a
            #       remote site, but not on the submit host. Currently allowed
            #       for the OSG site, versioning will not work as planned if
            #       we can't see the executable (can we perhaps run versioning
            #       including singularity??)

            # Check that executables at file urls
            #  on the local site exist
            if os.path.isfile(exe_url.path) is False:
                raise TypeError("Failed to find %s executable "
                                "at %s on site %s" % (name, exe_path,
                                exe_site))
        elif exe_url.scheme == 'singularity':
            # Will use an executable within a singularity container. Don't
            # need to do anything here, as I cannot easily check it exists.
            exe_path = exe_url.path
        else:
            # Could be http, https, etc. so it needs fetching if run now
            self.needs_fetching = True
            if self.needs_fetching and not self.installed:
                err_msg = "Non-file path URLs cannot be used unless the "
                err_msg += "executable is a bundled standalone executable. "
                err_msg += "If this is the case, then add the "
                err_msg += "pycbc.installed=True property."
                raise ValueError(err_msg)

        if self.installed:
            # Is installed, so copy from local site, like other inputs
            self.exe_pfns['local'] = exe_path
        else:
            # We must rely on the executables, and accompanying libraries,
            # being directly accessible on the execution site.
            # CVMFS is perfect for this! As is singularity.
            self.exe_pfns[exe_site] = exe_path
        logging.debug(
            "Using %s executable at %s on site %s",
            name,
            exe_url.path,
            exe_site
        )

        # FIXME: This hasn't yet been ported to pegasus5 and won't work.
        #        Pegasus describes two ways to work with containers, and I need
        #        to figure out which is most appropriate and use that.
        # Determine if this executables should be run in a container
        try:
            self.container_type = cp.get('pegasus_profile-%s' % name,
                                         'container|type')
        except:
            pass

        if self.container_type is not None:
            # FIXME: Move the actual container setup into pegasus_workflow
            self.container_img = cp.get('pegasus_profile-%s' % name,
                                        'container|image')
            try:
                self.container_site = cp.get('pegasus_profile-%s' % name,
                                             'container|image_site')
            except:
                self.container_site = 'local'

            try:
                self.container_mount = cp.get('pegasus_profile-%s' % name,
                                             'container|mount').split(',')
            except:
                self.container_mount = None


            self.container_cls = Pegasus.api.Container("{}-container".format(
                                                    name),
                                                    self.container_type,
                                                    self.container_img,
                                                    imagesite=self.container_site,
                                                    mount=self.container_mount)

            super(Executable, self).__init__(self.pegasus_name,
                                             installed=self.installed,
                                             container=self.container_cls)

        else:
            super(Executable, self).__init__(self.pegasus_name,
                                             installed=self.installed)

        if hasattr(self, "group_jobs"):
            self.add_profile('pegasus', 'clusters.size', self.group_jobs)

        # This sets up the sub-directory to use in the submit directory
        if set_submit_subdir:
            self.add_profile('pegasus', 'relative.submit.dir',
                             self.pegasus_name)

        # Set configurations from the config file, these should override all
        # other settings
        self._set_pegasus_profile_options()

        self.execution_site = exe_site
        self.executable_url = exe_path

    @property
    def ifo(self):
        """Return the ifo.

        If only one ifo in the ifo list this will be that ifo. Otherwise an
        error is raised.
        """
        if self.ifo_list and len(self.ifo_list) == 1:
            return self.ifo_list[0]
        else:
            errMsg = "self.ifoList must contain only one ifo to access the "
            errMsg += "ifo property. %s." %(str(self.ifo_list),)
            raise TypeError(errMsg)

    def get_transformation(self):
        if self.execution_site in self.transformations:
            return self.transformations[self.execution_site]
        else:
            self.create_transformation(self.execution_site,
                                       self.executable_url)
            return self.get_transformation()

    def add_ini_profile(self, cp, sec):
        """Add profile from configuration file.

        Parameters
        -----------
        cp : ConfigParser object
            The ConfigParser object holding the workflow configuration settings
        sec : string
            The section containing options for this job.
        """
        for opt in cp.options(sec):
            namespace = opt.split('|')[0]
            if namespace == 'pycbc' or namespace == 'container':
                continue

            value = cp.get(sec, opt).strip()
            key = opt.split('|')[1]
            self.add_profile(namespace, key, value)

    def _add_ini_opts(self, cp, sec, ignore_existing=False):
        """Add job-specific options from configuration file.

        Parameters
        -----------
        cp : ConfigParser object
            The ConfigParser object holding the workflow configuration
            settings
        sec : string
            The section containing options for this job.
        """
        for opt in cp.options(sec):
            if opt in self.all_added_options:
                if ignore_existing:
                    continue
                else:
                    raise ValueError("Option %s has already been added" % opt)
            self.all_added_options.add(opt)
            value = cp.get(sec, opt).strip()
            opt = f'--{opt}'
            if opt in self.file_input_options:
                # This now expects the option to be a file
                # Check if we have a list of files
                values = [path for path in value.split(' ') if path]

                self.common_raw_options.append(opt)
                self.common_raw_options.append(' ')

                # Get LFN and PFN
                for path in values:
                    # Here I decide if the path is URL or
                    # IFO:/path/to/file or IFO:url://path/to/file
                    # That's somewhat tricksy as we used : as delimiter
                    split_path = path.split(':', 1)
                    if len(split_path) == 1:
                        ifo = None
                        path = path
                    else:
                        # Have I split a URL or not?
                        if split_path[1].startswith('//'):
                            # URL
                            ifo = None
                            path = path
                        else:
                            #IFO:path or IFO:URL
                            ifo = split_path[0]
                            path = split_path[1]

                    # If the file exists make sure to use the
                    # fill path as a file:// URL
                    if os.path.isfile(path):
                        curr_pfn = urljoin('file:',
                                           pathname2url(os.path.abspath(path)))
                    else:
                        curr_pfn = path
                    curr_file = resolve_url_to_file(curr_pfn)
                    self.common_input_files.append(curr_file)
                    if ifo:
                        self.common_raw_options.append(ifo + ':')
                        self.common_raw_options.append(curr_file.dax_repr)
                    else:
                        self.common_raw_options.append(curr_file.dax_repr)
                    self.common_raw_options.append(' ')
                    
            elif opt in self.time_dependent_options:
                # There is a possibility of time-dependent, file options.
                # For now we will avoid supporting that complication unless
                # it is needed. This would require resolving the file first
                # in this function, and then dealing with the time-dependent
                # stuff later.
                self.unresolved_td_options[opt] = value
            else:
                # This option comes from the config file(s)
                self.common_options += [opt, value]
                
    def add_opt(self, opt, value=None):
        """Add option to job.

        Parameters
        -----------
        opt : string
            Name of option (e.g. --output-file-format)
        value : string, (default=None)
            The value for the option (no value if set to None).
        """
        if value is None:
            self.common_options += [opt]
        else:
            self.common_options += [opt, value]
    
    def get_opt(self, opt):
        """Get value of option from configuration file

        Parameters
        -----------
        opt : string
            Name of option (e.g. output-file-format)

        Returns
        --------
        value : string
            The value for the option. Returns None if option not present.
        """
        for sec in self.sections:
            try:
                key = self.cp.get(sec, opt)
                if key:
                    return key
            except ConfigParser.NoOptionError:
                pass

        return None

    def has_opt(self, opt):
        """Check if option is present in configuration file

        Parameters
        -----------
        opt : string
            Name of option (e.g. output-file-format)
        """
        for sec in self.sections:
            val = self.cp.has_option(sec, opt)
            if val:
                return val

        return False

    def create_node(self, **kwargs):
        """Default node constructor.

        This is usually overridden by subclasses of Executable.
        """
        return Node(self, **kwargs)

    def update_current_retention_level(self, value):
        """Set a new value for the current retention level.

        This updates the value of self.retain_files for an updated value of the
        retention level.

        Parameters
        -----------
        value : int
            The new value to use for the retention level.
        """
        # Determine the level at which output files should be kept
        self.current_retention_level = value
        try:
            global_retention_level = \
                self.cp.get_opt_tags("workflow", "file-retention-level",
                                   self.tags+[self.name])
        except ConfigParser.Error:
            msg="Cannot find file-retention-level in [workflow] section "
            msg+="of the configuration file. Setting a default value of "
            msg+="retain all files."
            logging.warn(msg)
            self.retain_files = True
            self.global_retention_threshold = 1
            self.cp.set("workflow", "file-retention-level", "all_files")
        else:
            # FIXME: Are these names suitably descriptive?
            retention_choices = {
                                 'all_files' : 1,
                                 'all_triggers' : 2,
                                 'merged_triggers' : 3,
                                 'results' : 4
                                }
            try:
                self.global_retention_threshold = \
                      retention_choices[global_retention_level]
            except KeyError:
                err_msg = "Cannot recognize the file-retention-level in the "
                err_msg += "[workflow] section of the ini file. "
                err_msg += "Got : {0}.".format(global_retention_level)
                err_msg += "Valid options are: 'all_files', 'all_triggers',"
                err_msg += "'merged_triggers' or 'results' "
                raise ValueError(err_msg)
            if self.current_retention_level == 5:
                self.retain_files = True
                if type(self).__name__ in Executable._warned_classes_list:
                    pass
                else:
                    warn_msg = "Attribute current_retention_level has not "
                    warn_msg += "been set in class {0}. ".format(type(self))
                    warn_msg += "This value should be set explicitly. "
                    warn_msg += "All output from this class will be stored."
                    logging.warn(warn_msg)
                    Executable._warned_classes_list.append(type(self).__name__)
            elif self.global_retention_threshold > self.current_retention_level:
                self.retain_files = False
            else:
                self.retain_files = True

    def update_current_tags(self, tags):
        """Set a new set of tags for this executable.

        Update the set of tags that this job will use. This updated default
        file naming and shared options. It will *not* update the pegasus
        profile, which belong to the executable and cannot be different for
        different nodes.

        Parameters
        -----------
        tags : list
            The new list of tags to consider.
        """
        if tags is None:
            tags = []
        if '' in tags:
            logging.warn('DO NOT GIVE ME EMPTY TAGS (in %s)', self.name)
            tags.remove('')
        tags = [tag.upper() for tag in tags]
        self.tags = tags

        if len(tags) > 6:
            warn_msg = "This job has way too many tags. "
            warn_msg += "Current tags are {}. ".format(' '.join(tags))
            warn_msg += "Current executable {}.".format(self.name)
            logging.info(warn_msg)

        if len(tags) != 0:
            self.tagged_name = "{0}-{1}".format(self.name, '_'.join(tags))
        else:
            self.tagged_name = self.name
        if self.ifo_string is not None:
            self.tagged_name = "{0}-{1}".format(self.tagged_name,
                                                self.ifo_string)

        # Determine the sections from the ini file that will configure
        # this executable
        sections = [self.name]
        if self.ifo_list is not None:
            if len(self.ifo_list) > 1:
                sec_tags = tags + self.ifo_list + [self.ifo_string]
            else:
                sec_tags = tags + self.ifo_list
        else:
            sec_tags = tags
        for sec_len in range(1, len(sec_tags)+1):
            for tag_permutation in permutations(sec_tags, sec_len):
                joined_name = '-'.join(tag_permutation)
                section = '{0}-{1}'.format(self.name, joined_name.lower())
                if self.cp.has_section(section):
                    sections.append(section)
                    
        # specfically for hierarchical search: This can be written in a better way 
        if tags[0] == str('FULL_DATA'):
            sections.append('inspiral-fulldata')
            
        elif tags[0] == str('INJ1'):
            sections.append('inspiral-injset1')  
           
        elif tags[0] == str('INJ2'):
            sections.append('inspiral-injset2')
            
        elif tags[0] == str('INJ3'):
            sections.append('inspiral-injset3')
        
        self.sections = sections
        
        # Do some basic sanity checking on the options
        for sec1, sec2 in combinations(sections, 2):
            self.cp.check_duplicate_options(sec1, sec2, raise_error=True)

        # collect the options and profile information
        # from the ini file section(s)
        self.all_added_options = set()
        self.common_options = []
        self.common_raw_options = []
        self.unresolved_td_options = {}
        self.common_input_files = []
        
        for sec in sections:
            if self.cp.has_section(sec):
                self._add_ini_opts(self.cp, sec)
            else:
                warn_string = "warning: config file is missing section "
                warn_string += "[{0}]".format(sec)
                logging.warn(warn_string)

        # get uppermost section
        if self.cp.has_section(f'{self.name}-defaultvalues'):
            self._add_ini_opts(self.cp, f'{self.name}-defaultvalues',
                               ignore_existing=True)
            
    def update_output_directory(self, out_dir=None):
        """Update the default output directory for output files.

        Parameters
        -----------
        out_dir : string (optional, default=None)
            If provided use this as the output directory. Else choose this
            automatically from the tags.
        """
        # Determine the output directory
        if out_dir is not None:
            self.out_dir = out_dir
        elif len(self.tags) == 0:
            self.out_dir = self.name + "_output"
        else:
            self.out_dir = self.tagged_name + "_output"

        if not os.path.isabs(self.out_dir):
            self.out_dir = os.path.join(os.getcwd(), self.out_dir)

        # Make output directory if not there
        if not os.path.isdir(self.out_dir):
            make_analysis_dir(self.out_dir)

    def _set_pegasus_profile_options(self):
        """Set the pegasus-profile settings for this Executable.

        These are a property of the Executable and not of nodes that it will
        spawn. Therefore it *cannot* be updated without also changing values
        for nodes that might already have been created. Therefore this is
        only called once in __init__. Second calls to this will fail.
        """
        # Executable- and tag-specific profile information
        for sec in self.sections:
            if self.cp.has_section('pegasus_profile-{0}'.format(sec)):
                self.add_ini_profile(self.cp,
                                     'pegasus_profile-{0}'.format(sec))
                
                
def resolve_url_to_file(curr_pfn, attrs=None):
    """
    Resolves a PFN into a workflow.File object.

    This function will resolve a PFN to a workflow.File object. If a File
    object already exists for that PFN that will be returned, otherwise a new
    object is returned. We will implement default site schemes here as needed,
    for example cvfms paths will be added to the osg and nonfsio sites in
    addition to local. If the LFN is a duplicate of an existing one, but with a
    different PFN an AssertionError is raised. The attrs keyword-argument can
    be used to specify attributes of a file. All files have 4 possible
    attributes. A list of ifos, an identifying string - usually used to give
    the name of the executable that created the file, a segmentlist over which
    the file is valid and tags specifying particular details about those files.
    If attrs['ifos'] is set it will be used as the ifos, otherwise this will
    default to ['H1', 'K1', 'L1', 'V1']. If attrs['exe_name'] is given this
    will replace the "exe_name" sent to File.__init__ otherwise 'INPUT' will
    be given. segs will default to [[1,2000000000]] unless overridden with
    attrs['segs']. tags will default to an empty list unless overriden
    with attrs['tag']. If attrs is None it will be ignored and all defaults
    will be used. It is emphasized that these attributes are for the most part
    not important with input files. Exceptions include things like input
    template banks, where ifos and valid times will be checked in the workflow
    and used in the naming of child job output files.
    """
    cvmfsstr1 = 'file:///cvmfs/'
    cvmfsstr2 = 'file://localhost/cvmfs/'
    osdfstr1 = 'osdf:///'  # Technically this isn't CVMFS, but same handling!
    cvmfsstrs = (cvmfsstr1, cvmfsstr2, osdfstr1)

    # Get LFN
    urlp = urllib.parse.urlparse(curr_pfn)
    curr_lfn = os.path.basename(urlp.path)
    # Does this already exist as a File?
    if curr_lfn in file_input_from_config_dict.keys():
        file_pfn = file_input_from_config_dict[curr_lfn][2]
        # If the PFNs are different, but LFNs are the same then fail.
        assert(file_pfn == curr_pfn)
        curr_file = file_input_from_config_dict[curr_lfn][1]
    else:
        # Use resolve_url to download file/symlink as appropriate
        local_file_path = resolve_url(curr_pfn)
        # Create File object with default local path
        curr_file = File.from_path(local_file_path, attrs=attrs)
        if curr_pfn.startswith(cvmfsstrs):
            # Add PFNs for nonlocal sites for special cases (e.g. CVMFS).
            # This block could be extended as needed
            curr_file.add_pfn(curr_pfn, site='all')
        else:
            pfn_local = urljoin('file:', pathname2url(local_file_path))
            curr_file.add_pfn(pfn_local, 'local')
        # Store the file to avoid later duplication
        tuple_val = (local_file_path, curr_file, curr_pfn)
        file_input_from_config_dict[curr_lfn] = tuple_val
    return curr_file

class File(pegasus_workflow.File):
    '''
    This class holds the details of an individual output file
    This file(s) may be pre-supplied, generated from within the workflow
    command line script, or generated within the workflow. The important stuff
    is:

    * The ifo that the File is valid for
    * The time span that the OutFile is valid for
    * A short description of what the file is
    * The extension that the file should have
    * The url where the file should be located

    An example of initiating this class:

    >> c = File("H1", "INSPIRAL_S6LOWMASS", segments.segment(815901601, 815902001), file_url="file://localhost/home/spxiwh/H1-INSPIRAL_S6LOWMASS-815901601-400.xml.gz" )

    another where the file url is generated from the inputs:

    >> c = File("H1", "INSPIRAL_S6LOWMASS", segments.segment(815901601, 815902001), directory="/home/spxiwh", extension="xml.gz" )

    '''
    def __init__(self, ifos, exe_name, segs, file_url=None,
                 extension=None, directory=None, tags=None,
                 store_file=True, use_tmp_subdirs=False):
        """
        Create a File instance

        Parameters
        ----------
        ifos : string or list
            The ifo(s) that the File is valid for. If the file is
            independently valid for multiple ifos it can be provided as a list.
            Ie. ['H1',L1','V1'], if the file is only valid for the combination
            of ifos (for e.g. ligolw_thinca output) then this can be supplied
            as, for e.g. "H1L1V1".
        exe_name: string
            A short description of the executable description, tagging
            only the program that ran this job.
        segs : ligo.segments.segment or ligo.segments.segmentlist
            The time span that the OutFile is valid for. Note that this is
            *not* the same as the data that the job that made the file reads in.
            Lalapps_inspiral jobs do not analyse the first an last 72s of the
            data that is read, and are therefore not valid at those times. If
            the time is not continuous a segmentlist can be supplied.
        file_url : url (optional, default=None)
            If this is *not* supplied, extension and directory must be given.
            If specified this explicitly points to the url of the file, or the
            url where the file will be generated when made in the workflow.
        extension : string (optional, default=None)
            Either supply this *and* directory *or* supply only file_url.
            If given this gives the extension at the end of the file name. The
            full file name will be inferred from the other arguments
            following the workflow standard.
        directory : string (optional, default=None)
            Either supply this *and* extension *or* supply only file_url.
            If given this gives the directory in which the file exists, or will
            exists. The file name will be inferred from the other arguments
            following the workflow standard.
        tags : list of strings (optional, default=None)
            This is a list of descriptors describing what this file is. For
            e.g. this might be ["BNSINJECTIONS" ,"LOWMASS","CAT_2_VETO"].
            These are used in file naming.
        """
        self.metadata = {}

        # Set the science metadata on the file
        if isinstance(ifos, str):
            self.ifo_list = [ifos]
        else:
            self.ifo_list = ifos

        if self.ifo_list is not None:
            self.ifo_string = ''.join(self.ifo_list)
        else:
            self.ifo_string = 'file'

        self.description = exe_name

        if isinstance(segs, segments.segment):
            self.segment_list = segments.segmentlist([segs])
        elif isinstance(segs, (segments.segmentlist)):
            self.segment_list = segs
        else:
            err = "segs input must be either ligo.segments.segment or "
            err += "segments.segmentlist. Got %s." %(str(type(segs)),)
            raise ValueError(err)
        if tags is None:
            tags = []
        if '' in tags:
            logging.warn('DO NOT GIVE EMPTY TAGS (from %s)', exe_name)
            tags.remove('')
        self.tags = tags

        if len(self.tags):
            self.tag_str = '_'.join(tags)
            tagged_description = '_'.join([self.description] + tags)
        else:
            tagged_description = self.description

        # Follow the capitals-for-naming convention
        self.ifo_string = self.ifo_string.upper()
        self.tagged_description = tagged_description.upper()

        if not file_url:
            if not extension:
                raise TypeError("a file extension required if a file_url "
                                "is not provided")
            if not directory:
                raise TypeError("a directory is required if a file_url is "
                                "not provided")

            filename = self._filename(self.ifo_string, self.tagged_description,
                                      extension, self.segment_list.extent())
            path = os.path.join(directory, filename)
            if not os.path.isabs(path):
                path = os.path.join(os.getcwd(), path)
            file_url = urllib.parse.urlunparse(['file', 'localhost', path,
                                                None, None, None])

        if use_tmp_subdirs and len(self.segment_list):
            pegasus_lfn = str(int(self.segment_list.extent()[0]))[:-4]
            pegasus_lfn = pegasus_lfn + '/' + os.path.basename(file_url)
        else:
            pegasus_lfn = os.path.basename(file_url)
        super(File, self).__init__(pegasus_lfn)

        if store_file:
            self.storage_path = urllib.parse.urlsplit(file_url).path
        else:
            self.storage_path = None

    def __getstate__(self):
        """ Allow the workflow.File to be picklable. This disables the usage of
        the internal cache entry.
        """
        for i, seg in enumerate(self.segment_list):
            self.segment_list[i] = segments.segment(float(seg[0]), float(seg[1]))
        self.cache_entry = None
        safe_dict = copy.copy(self.__dict__)
        safe_dict['cache_entry'] = None
        return safe_dict

    # FIXME: This is a pegasus_workflow thing (don't think it's needed at all!)
    #        use the pegasus function directly (maybe not).
    def add_metadata(self, key, value):
        """ Add arbitrary metadata to this file """
        self.metadata[key] = value

    @property
    def ifo(self):
        """
        If only one ifo in the ifo_list this will be that ifo. Otherwise an
        error is raised.
        """
        if len(self.ifo_list) == 1:
            return self.ifo_list[0]
        else:
            err = "self.ifo_list must contain only one ifo to access the "
            err += "ifo property. %s." %(str(self.ifo_list),)
            raise TypeError(err)

    @property
    def segment(self):
        """
        If only one segment in the segmentlist this will be that segment.
        Otherwise an error is raised.
        """
        if len(self.segment_list) == 1:
            return self.segment_list[0]
        else:
            err = "self.segment_list must only contain one segment to access"
            err += " the segment property. %s." %(str(self.segment_list),)
            raise TypeError(err)

    @property
    def cache_entry(self):
        """
        Returns a CacheEntry instance for File.
        """
        if self.storage_path is None:
            raise ValueError('This file is temporary and so a lal '
                             'cache entry cannot be made')

        file_url = urllib.parse.urlunparse(['file', 'localhost',
                                            self.storage_path, None,
                                            None, None])
        cache_entry = lal.utils.CacheEntry(self.ifo_string,
                   self.tagged_description, self.segment_list.extent(), file_url)
        cache_entry.workflow_file = self
        return cache_entry

    def _filename(self, ifo, description, extension, segment):
        """
        Construct the standard output filename. Should only be used internally
        of the File class.
        """

        # Follow the frame convention of using integer filenames,
        # but stretching to cover partially covered seconds.
        start = int(segment[0])
        end = int(math.ceil(segment[1]))
        duration = str(end-start)
        start = str(start)
        
        if extension.startswith('.'):
            extension = extension[1:]
            return "%s-%s-%s-%s-%s" % (ifo, description.upper(), start,
                                   duration, extension)
        else:
            extension = extension
            return "%s-%s-%s-%s-%s" % (ifo, description.upper(), start,
                                   duration, extension)
    
    @classmethod
    def from_path(cls, path, attrs=None, **kwargs):
        """
        Create an output File object from path, with optional attributes.
        """
        if attrs is None:
            attrs = {}
        if attrs and 'ifos' in attrs:
            ifos = attrs['ifos']
        else:
            ifos = ['H1', 'K1', 'L1', 'V1']
        if attrs and 'exe_name' in attrs:
            exe_name = attrs['exe_name']
        else:
            exe_name = 'INPUT'
        if attrs and 'segs' in attrs:
            segs = attrs['segs']
        else:
            segs = segments.segment([1, 2000000000])
        if attrs and 'tags' in attrs:
            tags = attrs['tags']
        else:
            tags = []

        curr_file = cls(ifos, exe_name, segs, path, tags=tags, **kwargs)
        return curr_file

class Node(pegasus_workflow.Node):
    def __init__(self, executable, valid_seg=None):
        super(Node, self).__init__(executable.get_transformation())
        self.executable = executable
        self.executed = False
        self.set_category(executable.name)
        self.valid_seg = valid_seg

        self._options += self.executable.common_options
        self._raw_options += self.executable.common_raw_options
        for inp in self.executable.common_input_files:
            self.add_input(inp)

        if len(self.executable.time_dependent_options):
            # Resolving these options requires the concept of a valid time.
            # To keep backwards compatibility we will allow this to work if
            # valid_seg is not supplied and no option actually needs resolving.
            # It would be good to get this from the workflow's valid_seg if
            # not overriden. But the Node is not connected to the Workflow
            # until the dax starts to be written.
            self.resolve_td_options(self.executable.unresolved_td_options)

    def get_command_line(self):
        # FIXME: Put in pegasus_workflow??
        self._finalize()
        arglist = self._dax_node.arguments

        tmpargs = []
        for a in arglist:
            if not isinstance(a, File):
                tmpargs += a.split(' ')
            else:
                tmpargs.append(a)
        arglist = tmpargs

        arglist = [a for a in arglist if a != '']

        arglist = [a.storage_path if isinstance(a, File) else a for a in arglist]

        # This allows the pfn to be an http(s) URL, which will be
        # downloaded by resolve_url
        exe_path = urllib.parse.urlsplit(self.executable.get_pfn()).path

        return [exe_path] + arglist

    def new_output_file_opt(self, valid_seg, extension, option_name, tags=None,
                            store_file=None, use_tmp_subdirs=False):
        """
        This function will create a workflow.File object corresponding to the given
        information and then add that file as output of this node.

        Parameters
        -----------
        valid_seg : ligo.segments.segment
            The time span over which the job is valid for.
        extension : string
            The extension to be used at the end of the filename.
            E.g. '.xml' or '.sqlite'.
        option_name : string
            The option that is used when setting this job as output. For e.g.
            'output-name' or 'output-file', whatever is appropriate for the
            current executable.
        tags : list of strings, (optional, default=[])
            These tags will be added to the list of tags already associated with
            the job. They can be used to uniquely identify this output file.
        store_file : Boolean, (optional, default=True)
            This file is to be added to the output mapper and will be stored
            in the specified output location if True. If false file will be
            removed when no longer needed in the workflow.
        """
        if tags is None:
            tags = []
        # Changing this from set(tags) to enforce order. It might make sense
        # for all jobs to have file names with tags in the same order.
        all_tags = copy.deepcopy(self.executable.tags)
        for tag in tags:
            if tag not in all_tags:
                all_tags.append(tag)
                
        store_file = store_file if store_file is not None else self.executable.retain_files
        fil = File(self.executable.ifo_list, self.executable.name,
                   valid_seg, extension=extension, store_file=store_file,
                   directory=self.executable.out_dir, tags=all_tags,
                   use_tmp_subdirs=use_tmp_subdirs)
        self.add_output_opt(option_name, fil)
        return fil

    def add_multiifo_input_list_opt(self, opt, inputs):
        """ Add an option that determines a list of inputs from multiple
            detectors. Files will be supplied as --opt ifo1:input1 ifo2:input2
            .....
        """
        # NOTE: Here we have to use the raw arguments functionality as the
        #       file and ifo are not space separated.
        self.add_raw_arg(opt)
        self.add_raw_arg(' ')
        for infile in inputs:
            self.add_raw_arg(infile.ifo)
            self.add_raw_arg(':')
            self.add_raw_arg(infile.name)
            self.add_raw_arg(' ')
            self.add_input(infile)

    def add_multiifo_output_list_opt(self, opt, outputs):
        """ Add an option that determines a list of outputs from multiple
            detectors. Files will be supplied as --opt ifo1:input1 ifo2:input2
            .....
        """
        # NOTE: Here we have to use the raw arguments functionality as the
        #       file and ifo are not space separated.
        self.add_raw_arg(opt)
        self.add_raw_arg(' ')
        for outfile in outputs:
            self.add_raw_arg(outfile.ifo)
            self.add_raw_arg(':')
            self.add_raw_arg(outfile.name)
            self.add_raw_arg(' ')
            self.add_output(outfile)

    def new_multiifo_output_list_opt(self, opt, ifos, analysis_time, extension,
                                     tags=None, store_file=None,
                                     use_tmp_subdirs=False):
        """ Add an option that determines a list of outputs from multiple
            detectors. Files will be supplied as --opt ifo1:input1 ifo2:input2
            .....
            File names are created internally from the provided extension and
            analysis time.
        """
        if tags is None:
            tags = []
        all_tags = copy.deepcopy(self.executable.tags)
        for tag in tags:
            if tag not in all_tags:
                all_tags.append(tag)

        output_files = FileList([])
        store_file = store_file if store_file is not None \
                                              else self.executable.retain_files
        
        for ifo in ifos:
            curr_file = File(ifo, self.executable.name, analysis_time,
                             extension=extension, store_file=store_file,
                             directory=self.executable.out_dir, tags=all_tags,
                             use_tmp_subdirs=use_tmp_subdirs)
            output_files.append(curr_file)
        self.add_multiifo_output_list_opt(opt, output_files)

    def resolve_td_options(self, td_options):
        for opt in td_options:
            new_opt = resolve_td_option(td_options[opt], self.valid_seg)
            self._options += [opt, new_opt]

    @property
    def output_files(self):
        return FileList(self._outputs)

    @property
    def output_file(self):
        """
        If only one output file return it. Otherwise raise an exception.
        """
        out_files = self.output_files
        if len(out_files) != 1:
            err_msg = "output_file property is only valid if there is a single"
            err_msg += " output file. Here there are "
            err_msg += "%d output files." %(len(out_files))
            raise ValueError(err_msg)
        return out_files[0]
    

class FileList(list):
    '''
    This class holds a list of File objects. It inherits from the
    built-in list class, but also allows a number of features. ONLY
    pycbc.workflow.File instances should be within a FileList instance.
    '''
    entry_class = File

    def categorize_by_attr(self, attribute):
        '''
        Function to categorize a FileList by a File object
        attribute (eg. 'segment', 'ifo', 'description').

        Parameters
        -----------
        attribute : string
           File object attribute to categorize FileList

        Returns
        --------
        keys : list
           A list of values for an attribute
        groups : list
           A list of FileLists
        '''

        # need to sort FileList otherwise using groupby without sorting does
        # 'AAABBBCCDDAABB' -> ['AAA','BBB','CC','DD','AA','BB']
        # and using groupby with sorting does
        # 'AAABBBCCDDAABB' -> ['AAAAA','BBBBB','CC','DD']
        flist = sorted(self, key=attrgetter(attribute), reverse=True)

        # use groupby to create lists
        groups = []
        keys = []
        for k, g in groupby(flist, attrgetter(attribute)):
            groups.append(FileList(g))
            keys.append(k)

        return keys, groups

    def find_output(self, ifo, time):
        '''Returns one File most appropriate at the given time/time range.

        Return one File that covers the given time, or is most
        appropriate for the supplied time range.

        Parameters
        -----------
        ifo : string
           Name of the ifo (or ifos) that the file should be valid for.
        time : int/float/LIGOGPStime or tuple containing two values
           If int/float/LIGOGPStime (or similar may of specifying one time) is
           given, return the File corresponding to the time. This calls
           self.find_output_at_time(ifo,time).
           If a tuple of two values is given, return the File that is
           **most appropriate** for the time range given. This calls
           self.find_output_in_range

        Returns
        --------
        pycbc_file : pycbc.workflow.File instance
           The File that corresponds to the time or time range
        '''
        # Determine whether I have a specific time, or a range of times
        try:
            lenTime = len(time)
        except TypeError:
            # This is if I have a single time
            outFile = self.find_output_at_time(ifo,time)
        else:
            # This is if I have a range of times
            if lenTime == 2:
                outFile = self.find_output_in_range(ifo,time[0],time[1])
            # This is if I got a list that had more (or less) than 2 entries
            if len(time) != 2:
                raise TypeError("I do not understand the input variable time")
        return outFile

    def find_output_at_time(self, ifo, time):
        '''
        Return File that covers the given time.

        Parameters
        -----------
        ifo : string
           Name of the ifo (or ifos) that the File should correspond to
        time : int/float/LIGOGPStime
           Return the Files that covers the supplied time. If no
           File covers the time this will return None.

        Returns
        --------
        list of File classes
           The Files that corresponds to the time.
         '''
        # Get list of Files that overlap time, for given ifo
        outFiles = [i for i in self if ifo in i.ifo_list and time in i.segment_list]
        if len(outFiles) == 0:
            # No OutFile at this time
            return None
        elif len(outFiles) == 1:
            # 1 OutFile at this time (good!)
            return outFiles
        else:
            # Multiple output files. Currently this is valid, but we may want
            # to demand exclusivity later, or in certain cases. Hence the
            # separation.
            return outFiles

    def find_outputs_in_range(self, ifo, current_segment, useSplitLists=False):
        """
        Return the list of Files that is most appropriate for the supplied
        time range. That is, the Files whose coverage time has the
        largest overlap with the supplied time range.

        Parameters
        -----------
        ifo : string
           Name of the ifo (or ifos) that the File should correspond to
        current_segment : ligo.segments.segment
           The segment of time that files must intersect.

        Returns
        --------
        FileList class
           The list of Files that are most appropriate for the time range
        """
        currsegment_list = segments.segmentlist([current_segment])

        # Get all files overlapping the window
        overlap_files = self.find_all_output_in_range(ifo, current_segment,
                                                    useSplitLists=useSplitLists)

        # By how much do they overlap?
        overlap_windows = [abs(i.segment_list & currsegment_list) for i in overlap_files]

        if not overlap_windows:
            return []

        # Return the File with the biggest overlap
        # Note if two File have identical overlap, the first is used
        # to define the valid segment
        overlap_windows = numpy.array(overlap_windows, dtype = int)
        segmentLst = overlap_files[overlap_windows.argmax()].segment_list

        # Get all output files with the exact same segment definition
        output_files = [f for f in overlap_files if f.segment_list==segmentLst]
        return output_files

    def find_output_in_range(self, ifo, start, end):
        '''
        Return the File that is most appropriate for the supplied
        time range. That is, the File whose coverage time has the
        largest overlap with the supplied time range. If no Files
        overlap the supplied time window, will return None.

        Parameters
        -----------
        ifo : string
           Name of the ifo (or ifos) that the File should correspond to
        start : int/float/LIGOGPStime
           The start of the time range of interest.
        end : int/float/LIGOGPStime
           The end of the time range of interest

        Returns
        --------
        File class
           The File that is most appropriate for the time range
        '''
        currsegment_list = segments.segmentlist([segments.segment(start, end)])

        # First filter Files corresponding to ifo
        outFiles = [i for i in self if ifo in i.ifo_list]

        if len(outFiles) == 0:
            # No OutFiles correspond to that ifo
            return None
        # Filter OutFiles to those overlapping the given window
        currSeg = segments.segment([start,end])
        outFiles = [i for i in outFiles \
                                  if i.segment_list.intersects_segment(currSeg)]

        if len(outFiles) == 0:
            # No OutFile overlap that time period
            return None
        elif len(outFiles) == 1:
            # One OutFile overlaps that period
            return outFiles[0]
        else:
            overlap_windows = [abs(i.segment_list & currsegment_list) \
                                                        for i in outFiles]
            # Return the File with the biggest overlap
            # Note if two File have identical overlap, this will return
            # the first File in the list
            overlap_windows = numpy.array(overlap_windows, dtype = int)
            return outFiles[overlap_windows.argmax()]
        
    def find_all_output_in_range(self, ifo, currSeg, useSplitLists=False):
        """
        Return all files that overlap the specified segment.
        """
        if not useSplitLists:
            # Slower, but simpler method
            outFiles = [i for i in self if ifo in i.ifo_list]
            outFiles = [i for i in outFiles
                        if i.segment_list.intersects_segment(currSeg)]
        else:
            # Faster, but more complicated
            # Basically only check if a subset of files intersects_segment by
            # using a presorted list. Sorting only happens once.
            if not self._check_split_list_validity():
                # FIXME: DO NOT hard code this.
                self._temporal_split_list(100)
            startIdx = int((currSeg[0] - self._splitListsStart) /
                           self._splitListsStep)
            # Add some small rounding here
            endIdx = (currSeg[1] - self._splitListsStart) / self._splitListsStep
            endIdx = int(endIdx - 0.000001)

            outFiles = []
            for idx in range(startIdx, endIdx + 1):
                if idx < 0 or idx >= self._splitListsNum:
                    continue
                outFilesTemp = [i for i in self._splitLists[idx]
                                if ifo in i.ifo_list]
                outFiles.extend([i for i in outFilesTemp
                                 if i.segment_list.intersects_segment(currSeg)])
                # Remove duplicates
                outFiles = list(set(outFiles))

        return self.__class__(outFiles)
    
    def find_output_with_tag(self, tag):
        """
        Find all files who have tag in self.tags
        """
        # Enforce upper case
        tag = tag.upper()
        return FileList([i for i in self if tag in i.tags])

    def find_output_without_tag(self, tag):
        """
        Find all files who do not have tag in self.tags
        """
        # Enforce upper case
        tag = tag.upper()
        return FileList([i for i in self if tag not in i.tags])

    def find_output_with_ifo(self, ifo):
        """
        Find all files who have ifo = ifo
        """
        # Enforce upper case
        ifo = ifo.upper()
        return FileList([i for i in self if ifo in i.ifo_list])

    def get_times_covered_by_files(self):
        """
        Find the coalesced intersection of the segments of all files in the
        list.
        """
        times = segments.segmentlist([])
        for entry in self:
            times.extend(entry.segment_list)
        times.coalesce()
        return times

    def convert_to_lal_cache(self):
        """
        Return all files in this object as a glue.lal.Cache object
        """
        lal_cache = gluelal.Cache([])
        for entry in self:
            try:
                lal_cache.append(entry.cache_entry)
            except ValueError:
                pass
        return lal_cache

    def _temporal_split_list(self,numSubLists):
        """
        This internal function is used to speed the code up in cases where a
        number of operations are being made to determine if files overlap a
        specific time. Normally such operations are done on *all* entries with
        *every* call. However, if we predetermine which files are at which
        times, we can avoid testing *every* file every time.

        We therefore create numSubLists distinct and equal length time windows
        equally spaced from the first time entry in the list until the last.
        A list is made for each window and files are added to lists which they
        overlap.

        If the list changes it should be captured and these split lists become
        invalid. Currently the testing for this is pretty basic
        """
        # Assume segment lists are coalesced!
        startTime = float( min([i.segment_list[0][0] for i in self]))
        endTime = float( max([i.segment_list[-1][-1] for i in self]))
        step = (endTime - startTime) / float(numSubLists)

        # Set up storage
        self._splitLists = []
        for idx in range(numSubLists):
            self._splitLists.append(FileList([]))

        # Sort the files

        for currFile in self:
            segExtent = currFile.segment_list.extent()
            startIdx = (segExtent[0] - startTime) / step
            endIdx = (segExtent[1] - startTime) / step
            # Add some small rounding here
            startIdx = int(startIdx - 0.001)
            endIdx = int(endIdx + 0.001)

            if startIdx < 0:
                startIdx = 0
            if endIdx >= numSubLists:
                endIdx = numSubLists - 1

            for idx in range(startIdx, endIdx + 1):
                self._splitLists[idx].append(currFile)

        # Set information needed to detect changes and to be used elsewhere
        self._splitListsLength = len(self)
        self._splitListsNum = numSubLists
        self._splitListsStart = startTime
        self._splitListsEnd = endTime
        self._splitListsStep = step
        self._splitListsSet = True

    def _check_split_list_validity(self):
        """
        See _temporal_split_list above. This function checks if the current
        split lists are still valid.
        """
        # FIXME: Currently very primitive, but needs to be fast
        if not (hasattr(self,"_splitListsSet") and (self._splitListsSet)):
            return False
        elif len(self) != self._splitListsLength:
            return False
        else:
            return True

    @classmethod
    def load(cls, filename):
        """
        Load a FileList from a pickle file
        """
        f = open(filename, 'r')
        return pickle.load(f)

    def dump(self, filename):
        """
        Output this FileList to a pickle file
        """
        f = open(filename, 'w')
        pickle.dump(self, f)

    def to_file_object(self, name, out_dir):
        """Dump to a pickle file and return an File object reference

        Parameters
        ----------
        name : str
            An identifier of this file. Needs to be unique.
        out_dir : path
            path to place this file

        Returns
        -------
        file : AhopeFile
        """
        make_analysis_dir(out_dir)

        file_ref = File('ALL', name, self.get_times_covered_by_files(),
                             extension='.pkl', directory=out_dir)
        self.dump(file_ref.storage_path)
        return file_ref


def resolve_td_option(val_str, valid_seg):
    """
    Take an option which might be time-dependent and resolve it

    Some options might take different values depending on the GPS time. For
    example if you want opt_1 to take value_a if the time is between 10 and
    100, value_b if between 100 and 250, and value_c if between 250 and 500 you
    can supply:

    value_a[10:100],value_b[100:250],value_c[250:500].

    This function will parse that string (as opt) and return the value fully
    contained in valid_seg. If valid_seg is not full contained in one, and only
    one, of these options. The code will fail. If given a simple option like:

    value_a

    The function will just return value_a.
    """
    # Track if we've already found a matching option
    output = ''
    # Strip any whitespace, and split on comma
    curr_vals = val_str.replace(' ', '').strip().split(',')

    # Resolving the simple case is trivial and can be done immediately.
    if len(curr_vals) == 1 and '[' not in curr_vals[0]:
        return curr_vals[0]

    # Loop over all possible values
    for cval in curr_vals:
        start = int(valid_seg[0])
        end = int(valid_seg[1])
        # Extract limits for each case, and check overlap with valid_seg
        if '[' in cval:
            bopt = cval.split('[')[1].split(']')[0]
            start, end = bopt.split(':')
            cval = cval.replace('[' + bopt + ']', '')
        curr_seg = segments.segment(int(start), int(end))
        # The segments module is a bit weird so we need to check if the two
        # overlap using the following code. If valid_seg is fully within
        # curr_seg this will be true.
        if curr_seg.intersects(valid_seg) and \
                (curr_seg & valid_seg == valid_seg):
            if output:
                err_msg = "Time-dependent options must be disjoint."
                raise ValueError(err_msg)
            output = cval
    if not output:
        err_msg = "Could not resolve option {}".format(val_str)
        raise ValueError
    return output
