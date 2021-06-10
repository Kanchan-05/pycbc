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
import sys, os, stat, subprocess, logging, math, string
from six.moves import configparser as ConfigParser
from six.moves import urllib
from six.moves.urllib.request import pathname2url
from six.moves.urllib.parse import urljoin
#from six.moves import cPickle
#import copy
#import numpy, random
#from itertools import combinations, groupby, permutations
from itertools import combinations, permutations
#from operator import attrgetter
from six import string_types
#import lal
#import lal.utils
import Pegasus.DAX3
#from glue import lal as gluelal
from ligo import segments
#from glue.ligolw import table, lsctables, ligolw
#from glue.ligolw import utils as ligolw_utils
#from glue.ligolw.utils import segments as ligolw_segments
#from glue.ligolw.utils import process as ligolw_process
#from pycbc import makedir
#from pycbc.workflow.configuration import WorkflowConfigParser, resolve_url

from pycbc.workflow.configuration import resolve_url
from pycbc.workflow import pegasus_workflow
from pycbc.workflow.core import File

#class ContentHandler(ligolw.LIGOLWContentHandler):
#    pass

#lsctables.use_in(ContentHandler)

#REMOVE THESE FUNCTIONS  FOR PYTHON >= 2.7 ####################################
def check_output_error_and_retcode(*popenargs, **kwargs):
    """
    This function is used to obtain the stdout of a command. It is only used
    internally, recommend using the make_external_call command if you want
    to call external executables.
    """
    if 'stdout' in kwargs:
        raise ValueError('stdout argument not allowed, it will be overridden.')
    process = subprocess.Popen(stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               *popenargs, **kwargs)
    output, error = process.communicate()
    retcode = process.poll()
    return output, error, retcode

def check_output(*popenargs, **kwargs):
    output, _, _ = check_output_error_and_retcode(*popenargs, **kwargs)
    return output

###############################################################################

#def make_analysis_dir(path):
#    """
#    Make the analysis directory path, any parent directories that don't already
#    exist, and the 'logs' subdirectory of path.
#    """
#    if path is not None:
#        makedir(os.path.join(path, 'logs'))

def is_condor_exec(exe_path):
    """
    Determine if an executable is condor-compiled

    Parameters
    ----------
    exe_path : str
          The executable path

    Returns
    -------
    truth_value  : boolean
        Return True if the exe is condor compiled, False otherwise.
    """
    if str(check_output(['nm', '-a', exe_path])).find('condor') != -1:
        return True
    else:
        return False

file_input_from_config_dict = {}

# Executable is modified for h-search
class Executable(pegasus_workflow.Executable):
    # These are the file retention levels
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
    # file_input_options = ['--psd-file, '--bank-file'] (as an example)
    file_input_options = []
    
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
    def __init__(self, cp, name,
                 universe=None, ifos=None, out_dir=None, tags=None,
                 reuse_executable=True):
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
        
       
        if isinstance(ifos, string_types):
            self.ifo_list = [ifos]
        else:
            self.ifo_list = ifos
        if self.ifo_list is not None:
            self.ifo_list = sorted(self.ifo_list)
            self.ifo_string = ''.join(self.ifo_list)
        else:
            self.ifo_string = None
        self.cp = cp
        self.universe=universe
        self.container_cls = None
        self.container_type = None

        try:
            self.installed = cp.getboolean('pegasus_profile-%s' % name, 'pycbc|installed')
        except:
            self.installed = True

        self.name=name

        self.update_current_tags(tags)

        self.update_output_directory(out_dir=out_dir)

        # Determine the level at which output files should be kept
        self.update_current_retention_level(self.current_retention_level)

        # Should I reuse this executable?
        if reuse_executable:
            self.pegasus_name = self.name
        else:
            self.pegasus_name = self.tagged_name
            
        # Determine if this executables should be run in a container
        try:
            self.container_type = cp.get('pegasus_profile-%s' % name,
                                         'container|type')
        except:
            pass

        if self.container_type is not None:
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


            self.container_cls = Pegasus.DAX3.Container("{}-container".format(
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

        self._set_pegasus_profile_options()
        # Check that the executable actually exists locally or
        # looks like a URL, in which case trust Pegasus to be
        # able to fetch it.
        exe_path = cp.get('executables', name)
        self.needs_fetching = False

        exe_url = urllib.parse.urlparse(exe_path)

        # See if the user specified a list of sites for the executable
        try:
            exe_site_list = cp.get('pegasus_profile-%s' % name, 'pycbc|site')
        except:
            exe_site_list = 'local'

        for s in exe_site_list.split(','):
            exe_site = s.strip()

            if exe_url.scheme in ['', 'file']:
                if exe_site is 'local':
                    # Check that executables at file urls
                    #  on the local site exist
                    if os.path.isfile(exe_url.path) is False:
                        raise TypeError("Failed to find %s executable "
                                        "at %s on site %s" % (name, exe_path,
                                        exe_site))
            else:
                # Could be http, gsiftp, etc. so it needs fetching if run now
                self.needs_fetching = True

            self.add_pfn(exe_path, site=exe_site)
            logging.info("Using %s executable "
                         "at %s on site %s" % (name, exe_url.path, exe_site))

        # Determine the condor universe if we aren't given one
        if self.universe is None:
            if is_condor_exec(exe_path):
                self.universe = 'standard'
            else:
                self.universe = 'vanilla'

        if not self.universe == 'vanilla':
            logging.info("%s executable will run as %s universe"
                         % (name, self.universe))

        self.set_universe(self.universe)

        if hasattr(self, "group_jobs"):
            self.add_profile('pegasus', 'clusters.size', self.group_jobs)

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
            self.add_profile(namespace, key, value, force=True)

            # Remove if Pegasus can apply this hint in the TC
            if namespace == 'hints' and key == 'execution.site':
                self.execution_site = value

    def add_ini_opts(self, cp, sec):
        """Add job-specific options from configuration file.

        Parameters
        -----------
        cp : ConfigParser object
            The ConfigParser object holding the workflow configuration settings
        sec : string
            The section containing options for this job.
        """
        for opt in cp.options(sec):
            value = cp.get(sec, opt).strip()
            opt = '--%s' %(opt,)
            if opt in self.file_input_options:
                # This now expects the option to be a file
                # Check is we have a list of files
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
            logging.warn('DO NOT GIVE ME EMPTY TAGS')
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
       
        if tags[0] == str('FULL_DATA'):
            sections.append('inspiral-fulldata')
            
        elif tags[0] == str('BNSSTT2_INJ'):
            sections.append('inspiral-bnsinj')  
           
        elif tags[0] == str('BBHSEOBNRV4_INJ'):
            sections.append('inspiral-bbhinj')
            
        elif tags[0] == str('NSBHSEOBNRV4_INJ'):
            sections.append('inspiral-nsbhinj')
        
        self.sections = sections
        
        # Do some basic sanity checking on the options
        for sec1, sec2 in combinations(sections, 2):
            self.cp.check_duplicate_options(sec1, sec2, raise_error=True)
        
        
        # collect the options and profile information
        # from the ini file section(s)
        self.common_options = []
        self.common_raw_options = []
        self.unresolved_td_options = {}
        self.common_input_files = []
        for sec in sections:
            if self.cp.has_section(sec):
                self.add_ini_opts(self.cp, sec)
            else:
                warn_string = "warning: config file is missing section "
                warn_string += "[{0}]".format(sec)
                logging.warn(warn_string)

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
            self.out_dir = self.name
        else:
            self.out_dir = self.tagged_name
        if not os.path.isabs(self.out_dir):
            self.out_dir = os.path.join(os.getcwd(), self.out_dir)

    def _set_pegasus_profile_options(self):
        """Set the pegasus-profile settings for this Executable.

        These are a property of the Executable and not of nodes that it will
        spawn. Therefore it *cannot* be updated without also changing values
        for nodes that might already have been created. Therefore this is
        only called once in __init__. Second calls to this will fail.
        """
        # Add executable non-specific profile information
        if self.cp.has_section('pegasus_profile'):
            self.add_ini_profile(self.cp, 'pegasus_profile')

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
    cvmfsstrs = (cvmfsstr1, cvmfsstr2)

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
        # To do this we first need to check the attributes
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

        curr_file = File(ifos, exe_name, segs, local_file_path, tags=tags)
        pfn_local = urljoin('file:', pathname2url(local_file_path))
        curr_file.PFN(pfn_local, 'local')
        # Add other PFNs for nonlocal sites as needed.
        # This block could be extended as needed
        if curr_pfn.startswith(cvmfsstrs):
            curr_file.PFN(curr_pfn, site='osg')
            curr_file.PFN(curr_pfn, site='nonfsio')
            # Also register the CVMFS PFN with the local site. We want to
            # prefer this, and symlink from here, when possible.
            # However, I think we need a little more to avoid it symlinking
            # to this through an NFS mount.
            curr_file.PFN(curr_pfn, site='local')
        # Store the file to avoid later duplication
        tuple_val = (local_file_path, curr_file, curr_pfn)
        file_input_from_config_dict[curr_lfn] = tuple_val
    return curr_file

