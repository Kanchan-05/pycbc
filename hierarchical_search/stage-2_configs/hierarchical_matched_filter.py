from __future__ import division

import os, logging
import math
#import lal
#from math import radians
#import pycbc
from pycbc.workflow.core import make_analysis_dir, File, FileList, Node
from ligo import segments
#import Pegasus.DAX3 as dax
#from pycbc.workflow.core import Executable, File, FileList, Node
from pycbc.workflow.jobsetup import select_tmpltbank_class
from hierarchical_jobsetup import PyCBCHierarchicalInspiralExecutable

        
class JobSegmenter(object):
    """ This class is used when running sngl_ifo_job_setup to determine what times
    should be analysed be each job and what data is needed.
    """

    def __init__(self, data_lengths, valid_chunks, valid_lengths, curr_seg,
                 curr_exe_class, compatibility_mode = False):
        """ Initialize class. """
        self.exe_class = curr_exe_class
        self.curr_seg = curr_seg
        self.curr_seg_length = float(abs(curr_seg))

        self.data_length, self.valid_chunk, self.valid_length = \
                      self.pick_tile_size(self.curr_seg_length, data_lengths,
                                                  valid_chunks, valid_lengths)

        self.data_chunk = segments.segment([0, self.data_length])
        self.data_loss = self.data_length - abs(self.valid_chunk)

        if self.data_loss < 0:
            raise ValueError("pycbc.workflow.jobsetup needs fixing! Please contact a developer")

        if self.curr_seg_length < self.data_length:
            self.num_jobs = 0
            return

        # How many jobs do we need
        self.num_jobs = int( math.ceil( (self.curr_seg_length \
                                - self.data_loss) / float(self.valid_length) ))

        self.compatibility_mode = compatibility_mode
        if compatibility_mode and (self.valid_length != abs(self.valid_chunk)):
            errMsg = "In compatibility mode the template bank and matched-"
            errMsg += "filter jobs must read in the same amount of data."
            print(self.valid_length, self.valid_chunk)
            raise ValueError(errMsg)
        elif compatibility_mode and len(data_lengths) > 1:
            raise ValueError("Cannot enable compatibility mode tiling with "
                             "variable tile size")
        elif compatibility_mode:
            # What is the incremental shift between jobs
            self.job_time_shift = self.valid_length
        elif self.curr_seg_length == self.data_length:
            # If the segment length is identical to the data length then I
            # will have exactly 1 job!
            self.job_time_shift = 0
        else:
            # What is the incremental shift between jobs
            self.job_time_shift = (self.curr_seg_length - self.data_length) / \
                                   float(self.num_jobs - 1)

    def pick_tile_size(self, seg_size, data_lengths, valid_chunks, valid_lengths):
        """ Choose job tiles size based on science segment length """

        if len(valid_lengths) == 1:
            return data_lengths[0], valid_chunks[0], valid_lengths[0]
        else:
            # Pick the tile size that is closest to 1/3 of the science segment
            target_size = seg_size / 3
            pick, pick_diff = 0, abs(valid_lengths[0] - target_size)
            for i, size in enumerate(valid_lengths):
                if abs(size - target_size) < pick_diff:
                    pick, pick_diff  = i, abs(size - target_size)
            return data_lengths[pick], valid_chunks[pick], valid_lengths[pick]

    def get_valid_times_for_job(self, num_job, allow_overlap=True):
        """ Get the times for which this job is valid. """
        if self.compatibility_mode:
            return self.get_valid_times_for_job_legacy(num_job)
        else:
            return self.get_valid_times_for_job_workflow(num_job,
                                                   allow_overlap=allow_overlap)

    def get_valid_times_for_job_workflow(self, num_job, allow_overlap=True):
        """ Get the times for which the job num_job will be valid, using workflow's
        method.
        """
        # small factor of 0.0001 to avoid float round offs causing us to
        # miss a second at end of segments.
        shift_dur = self.curr_seg[0] + int(self.job_time_shift * num_job\
                                           + 0.0001)
        job_valid_seg = self.valid_chunk.shift(shift_dur)
        # If we need to recalculate the valid times to avoid overlap
        if not allow_overlap:
            data_per_job = (self.curr_seg_length - self.data_loss) / \
                           float(self.num_jobs)
            lower_boundary = num_job*data_per_job + \
                                 self.valid_chunk[0] + self.curr_seg[0]
            upper_boundary = data_per_job + lower_boundary
            # NOTE: Convert to int after calculating both boundaries
            # small factor of 0.0001 to avoid float round offs causing us to
            # miss a second at end of segments.
            lower_boundary = int(lower_boundary)
            upper_boundary = int(upper_boundary + 0.0001)
            if lower_boundary < job_valid_seg[0] or \
                    upper_boundary > job_valid_seg[1]:
                err_msg = ("Workflow is attempting to generate output "
                          "from a job at times where it is not valid.")
                raise ValueError(err_msg)
            job_valid_seg = segments.segment([lower_boundary,
                                              upper_boundary])
        return job_valid_seg

    def get_valid_times_for_job_legacy(self, num_job):
        """ Get the times for which the job num_job will be valid, using the method
        use in inspiral hipe.
        """
        # All of this should be integers, so no rounding factors needed.
        shift_dur = self.curr_seg[0] + int(self.job_time_shift * num_job)
        job_valid_seg = self.valid_chunk.shift(shift_dur)

        # If this is the last job, push the end back
        if num_job == (self.num_jobs - 1):
            dataPushBack = self.data_length - self.valid_chunk[1]
            job_valid_seg = segments.segment(job_valid_seg[0],
                                               self.curr_seg[1] - dataPushBack)

        return job_valid_seg

    def get_data_times_for_job(self, num_job):
        """ Get the data that this job will read in. """
        if self.compatibility_mode:
            job_data_seg =  self.get_data_times_for_job_legacy(num_job)
        else:
            job_data_seg =  self.get_data_times_for_job_workflow(num_job)

        # Sanity check that all data is used
        if num_job == 0:
            if job_data_seg[0] != self.curr_seg[0]:
                err= "Job is not using data from the start of the "
                err += "science segment. It should be using all data."
                raise ValueError(err)
        if num_job == (self.num_jobs - 1):
            if job_data_seg[1] != self.curr_seg[1]:
                err = "Job is not using data from the end of the "
                err += "science segment. It should be using all data."
                raise ValueError(err)

        if hasattr(self.exe_class, 'zero_pad_data_extend'):
            job_data_seg = self.exe_class.zero_pad_data_extend(job_data_seg,
                                                                 self.curr_seg)

        return job_data_seg
    
    def get_data_times_for_job_workflow(self, num_job):
        """ Get the data that this job will need to read in. """
        # small factor of 0.0001 to avoid float round offs causing us to
        # miss a second at end of segments.
        shift_dur = self.curr_seg[0] + int(self.job_time_shift * num_job\
                                           + 0.0001)
        job_data_seg = self.data_chunk.shift(shift_dur)
        return job_data_seg

    def get_data_times_for_job_legacy(self, num_job):
        """ Get the data that this job will need to read in. """
        # Should all be integers, so no rounding needed
        shift_dur = self.curr_seg[0] + int(self.job_time_shift * num_job)
        job_data_seg = self.data_chunk.shift(shift_dur)

        # If this is the last job, push the end back
        if num_job == (self.num_jobs - 1):
            dataPushBack = job_data_seg[1] - self.curr_seg[1]
            assert dataPushBack >= 0
            job_data_seg = segments.segment(job_data_seg[0] - dataPushBack,
                                                 self.curr_seg[1])
            assert (abs(job_data_seg) == self.data_length)
        return job_data_seg
    
def identify_needed_data(curr_exe_job, link_job_instance=None):
    """ This function will identify the length of data that a specific executable
    needs to analyse and what part of that data is valid (ie. inspiral doesn't
    analyse the first or last 64+8s of data it reads in).

    In addition you can supply a second job instance to "link" to, which will
    ensure that the two jobs will have a one-to-one correspondence (ie. one
    template bank per one matched-filter job) and the corresponding jobs will
    be "valid" at the same times.

    Parameters
    -----------
    curr_exe_job : Job
        An instance of the Job class that has a get_valid times method.
    link_job_instance : Job instance (optional),
        Coordinate the valid times with another executable.

    Returns
    --------
    dataLength : float
        The amount of data (in seconds) that each instance of the job must read
        in.
    valid_chunk : glue.segment.segment
        The times within dataLength for which that jobs output **can** be
        valid (ie. for inspiral this is (72, dataLength-72) as, for a standard
        setup the inspiral job cannot look for triggers in the first 72 or
        last 72 seconds of data read in.)
    valid_length : float
        The maximum length of data each job can be valid for. If not using
        link_job_instance this is abs(valid_segment), but can be smaller than
        that if the linked job only analyses a small amount of data (for e.g.).
    """
    # Set up the condorJob class for the current executable
    data_lengths, valid_chunks = curr_exe_job.get_valid_times()

    # Begin by getting analysis start and end, and start and end of time
    # that the output file is valid for
    valid_lengths = [abs(valid_chunk) for valid_chunk in valid_chunks]

    if link_job_instance:
        # FIXME: Should we remove this, after testing is complete??
        # EURGHH! What we are trying to do here is, if this option is given,
        # line up the template bank and inspiral jobs so that there is one
        # template bank for each inspiral job. This is proving a little messy
        # and probably still isn't perfect.

        # What data does the linked exe use?
        link_data_length, link_valid_chunk = link_job_instance.get_valid_times()
        if len(link_data_length) > 1 or len(valid_lengths) > 1:
            raise ValueError('Linking job instances for tiling is not supported'
                             ' between jobs that allow variable tile size')

        # What data is lost at start of both jobs? Take the maximum.
        start_data_loss = max(valid_chunks[0][0], link_valid_chunk[0][0])
        # What data is lost at end of both jobs? Take the maximum.
        end_data_loss = max(data_lengths[0] - valid_chunks[0][1],\
                            link_data_length[0] - link_valid_chunk[0][1])
        # Calculate valid_segments for both jobs based on the combined data
        # loss.

        valid_chunks[0] = segments.segment(start_data_loss, \
                                       data_lengths[0] - end_data_loss)
        link_valid_chunk = segments.segment(start_data_loss, \
                                       link_data_length[0] - end_data_loss)

        # The maximum valid length should be the minimum of the two
        link_valid_length = abs(link_valid_chunk)

        # Which one is now longer? Use this is valid_length
        if link_valid_length < valid_lengths[0]:
            valid_lengths[0] = link_valid_length

    return data_lengths, valid_chunks, valid_lengths  

def sngl_ifo_job_setup(workflow, ifo, out_files, curr_exe_job, science_segs,
                       datafind_outs, parents=None,
                       link_job_instance=None, allow_overlap=True,
                       compatibility_mode=True):
    """ This function sets up a set of single ifo jobs. A basic overview of how this
    works is as follows:

    * (1) Identify the length of data that each job needs to read in, and what
      part of that data the job is valid for.
    * START LOOPING OVER SCIENCE SEGMENTS
    * (2) Identify how many jobs are needed (if any) to cover the given science
      segment and the time shift between jobs. If no jobs continue.
    * START LOOPING OVER JOBS
    * (3) Identify the time that the given job should produce valid output (ie.
      inspiral triggers) over.
    * (4) Identify the data range that the job will need to read in to produce
      the aforementioned valid output.
    * (5) Identify all parents/inputs of the job.
    * (6) Add the job to the workflow
    * END LOOPING OVER JOBS
    * END LOOPING OVER SCIENCE SEGMENTS

    Parameters
    -----------
    workflow: pycbc.workflow.core.Workflow
        An instance of the Workflow class that manages the constructed workflow.
    ifo : string
        The name of the ifo to set up the jobs for
    out_files : pycbc.workflow.core.FileList
        The FileList containing the list of jobs. Jobs will be appended
        to this list, and it does not need to be empty when supplied.
    curr_exe_job : Job
        An instanced of the Job class that has a get_valid times method.
    science_segs : ligo.segments.segmentlist
        The list of times that the jobs should cover
    datafind_outs : pycbc.workflow.core.FileList
        The file list containing the datafind files.
    parents : pycbc.workflow.core.FileList (optional, kwarg, default=None)
        The FileList containing the list of jobs that are parents to
        the one being set up.
    link_job_instance : Job instance (optional),
        Coordinate the valid times with another Executable.
    allow_overlap : boolean (optional, kwarg, default = True)
        If this is set the times that jobs are valid for will be allowed to
        overlap. This may be desired for template banks which may have some
        overlap in the times they cover. This may not be desired for inspiral
        jobs, where you probably want triggers recorded by jobs to not overlap
        at all.
    compatibility_mode : boolean (optional,  kwarg, default = False)
        If given the jobs will be tiled in the same method as used in inspiral
        hipe. This requires that link_job_instance is also given. If not given
        workflow's methods are used.

    Returns
    --------
    out_files : pycbc.workflow.core.FileList
        A list of the files that will be generated by this step in the
        workflow.
    """
    if compatibility_mode and not link_job_instance:
        errMsg = "Compability mode requires a link_job_instance."
        raise ValueError(errMsg)

    ########### (1) ############
    # Get the times that can be analysed and needed data lengths
    data_length, valid_chunk, valid_length = identify_needed_data(curr_exe_job,
                                           link_job_instance=link_job_instance)

    # Loop over science segments and set up jobs
    for curr_seg in science_segs:
        ########### (2) ############
        # Initialize the class that identifies how many jobs are needed and the
        # shift between them.
        segmenter = JobSegmenter(data_length, valid_chunk, valid_length,
                                 curr_seg, curr_exe_job,
                                 compatibility_mode=compatibility_mode)

        for job_num in range(segmenter.num_jobs):
            ############## (3) #############
            # Figure out over what times this job will be valid for
            job_valid_seg = segmenter.get_valid_times_for_job(job_num,
                                                   allow_overlap=allow_overlap)

            ############## (4) #############
            # Get the data that this job should read in
            job_data_seg = segmenter.get_data_times_for_job(job_num)

            ############# (5) ############
            # Identify parents/inputs to the job
            if parents:
                # Find the set of files with the best overlap
                curr_parent = parents.find_outputs_in_range(ifo, job_valid_seg,
                                                            useSplitLists=True)
                if not curr_parent:
                    err_string = ("No parent jobs found overlapping %d to %d."
                                  %(job_valid_seg[0], job_valid_seg[1]))
                    err_string += "\nThis is a bad error! Contact a developer."
                    raise ValueError(err_string)
            else:
                curr_parent = [None]

            curr_dfouts = None
            if datafind_outs:
                curr_dfouts = datafind_outs.find_all_output_in_range(ifo,
                                              job_data_seg, useSplitLists=True)
                if not curr_dfouts:
                    err_str = ("No datafind jobs found overlapping %d to %d."
                                %(job_data_seg[0],job_data_seg[1]))
                    err_str += "\nThis shouldn't happen. Contact a developer."
                    raise ValueError(err_str)


            ############## (6) #############
            # Make node and add to workflow

            # Note if I have more than one curr_parent I need to make more than
            # one job. If there are no curr_parents it is set to [None] and I
            # make a single job. This catches the case of a split template bank
            # where I run a number of jobs to cover a single range of time.

            # Sort parent jobs to ensure predictable order
            sorted_parents = sorted(curr_parent,
                                    key=lambda fobj: fobj.tagged_description)
            for pnum, parent in enumerate(sorted_parents):
                if len(curr_parent) != 1:
                    tag = ["JOB%d" %(pnum,)]
                else:
                    tag = []
                # To ensure output file uniqueness I add a tag
                # We should generate unique names automatically, but it is a
                # pain until we can set the output names for all Executables
                node = curr_exe_job.create_node(job_data_seg, job_valid_seg,
                                                parent=parent,
                                                dfParents=curr_dfouts,
                                                tags=tag)
                workflow.add_node(node)
                curr_out_files = node.output_files
                # FIXME: Here we remove PSD files if they are coming through.
                #        This should be done in a better way. On to-do list.
                curr_out_files = [i for i in curr_out_files if 'PSD_FILE'\
                                                                 not in i.tags]
                out_files += curr_out_files

    return out_files
# no change till here

# Added the option for PyCBCHierarchicalInspiralExecutable class existing in hierarchical_jobsetp_class 
def select_matchedfilter_class(curr_exe):
    """ This function returns a class that is appropriate for setting up
    matched-filtering jobs within workflow.

    Parameters
    ----------
    curr_exe : string
        The name of the matched filter executable to be used.

    Returns
    --------
    exe_class : Sub-class of pycbc.workflow.core.Executable that holds utility
        functions appropriate for the given executable.  Instances of the class
        ('jobs') **must** have methods
        * job.create_node()
        and
        * job.get_valid_times(ifo, )
    """

    exe_to_class_map = {
        'pycbc_hierarchical_inspiral'  : PyCBCHierarchicalInspiralExecutable
        #'pycbc_inspiral'              : PyCBCInspiralExecutable,
        #'pycbc_inspiral_skymax'       : PyCBCInspiralExecutable,
        #'pycbc_multi_inspiral'        : PyCBCMultiInspiralExecutable,
    }
    try:
        return exe_to_class_map[curr_exe]
    except KeyError:
        # also conceivable to introduce a default class??
        raise NotImplementedError(
            "No job class exists for executable %s, exiting" % curr_exe)

        
# option in usual matched_filter.py
def setup_matchedfltr_dax_generated(workflow, science_segs, datafind_outs,
                                    tmplt_banks, output_dir,
                                    injection_file=None,
                                    tags=None, link_to_tmpltbank=False,
                                    compatibility_mode=False):
    '''
    Setup matched-filter jobs that are generated as part of the workflow.
    This
    module can support any matched-filter code that is similar in principle to
    lalapps_inspiral, but for new codes some additions are needed to define
    Executable and Job sub-classes (see jobutils.py).
    Parameters
    -----------
    workflow : pycbc.workflow.core.Workflow
        The Workflow instance that the coincidence jobs will be added to.
    science_segs : ifo-keyed dictionary of ligo.segments.segmentlist instances
        The list of times that are being analysed in this workflow.
    datafind_outs : pycbc.workflow.core.FileList
        An FileList of the datafind files that are needed to obtain the
        data used in the analysis.
    tmplt_banks : pycbc.workflow.core.FileList
        An FileList of the template bank files that will serve as input
        in this stage.
    output_dir : path
        The directory in which output will be stored.
    injection_file : pycbc.workflow.core.File, optional (default=None)
        If given the file containing the simulation file to be sent to these
        jobs on the command line. If not given no file will be sent.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['BNSINJECTIONS'] or
        ['NOINJECTIONANALYSIS']. This will be used in output names.
    link_to_tmpltbank : boolean, optional (default=True)
        If this option is given, the job valid_times will be altered so that there
        will be one inspiral file for every template bank and they will cover the
        same time span. Note that this option must also be given during template
        bank generation to be meaningful.
    Returns
    -------
    inspiral_outs : pycbc.workflow.core.FileList
        A list of output files written by this stage. This *will not* contain
        any intermediate products produced within this stage of the workflow.
        If you require access to any intermediate products produced at this
        stage you can call the various sub-functions directly.
    '''
    if tags is None:
        tags = []
    # Need to get the exe to figure out what sections are analysed, what is
    # discarded etc. This should *not* be hardcoded, so using a new executable
    # will require a bit of effort here ....

    cp = workflow.cp
    ifos = science_segs.keys()
    print 'here---->>',tags
    match_fltr_exe = os.path.basename(cp.get('executables','inspiral'))
    # Select the appropriate class
    exe_class = select_matchedfilter_class(match_fltr_exe)

    if link_to_tmpltbank:
        # Use this to ensure that inspiral and tmpltbank jobs overlap. This
        # means that there will be 1 inspiral job for every 1 tmpltbank and
        # the data read in by both will overlap as much as possible. (If you
        # ask the template bank jobs to use 2000s of data for PSD estimation
        # and the matched-filter jobs to use 4000s, you will end up with
        # twice as many matched-filter jobs that still use 4000s to estimate a
        # PSD but then only generate triggers in the 2000s of data that the
        # template bank jobs ran on.
        tmpltbank_exe = os.path.basename(cp.get('executables', 'tmpltbank'))
        link_exe_instance = select_tmpltbank_class(tmpltbank_exe)
    else:
        link_exe_instance = None

    # Set up class for holding the banks
    inspiral_outs = FileList([])

    # Matched-filtering is done independently for different ifos, but might not be!
    # If we want to use multi-detector matched-filtering or something similar to this
    # it would probably require a new module
    for ifo in ifos:
        logging.info("Setting up matched-filtering for %s." %(ifo))
        job_instance = exe_class(workflow.cp,'inspiral',ifo=ifo,
                                               out_dir=output_dir,
                                               injection_file=injection_file,
                                               tags=tags)
        if link_exe_instance:
            link_job_instance = link_exe_instance(cp, 'tmpltbank', ifo=ifo,
                                               out_dir=output_dir, tags=tags)
        else:
            link_job_instance = None

        # changed specific to hierarchcial search---> kanchan
        sngl_ifo_job_setup(workflow, ifo, inspiral_outs, job_instance,
                           science_segs[ifo], datafind_outs,
                           parents=tmplt_banks, allow_overlap=False,
                           link_job_instance=link_job_instance,
                           compatibility_mode=compatibility_mode)
    return inspiral_outs

# we haven't added the option of setup_matchedfltr_dax_generated_multi as per matched_filter.py

# function in usual pycbc matched_filter.py
def setup_matchedfltr_workflow(workflow, science_segs, datafind_outs,
                               tmplt_banks, output_dir=None,
                               injection_file=None, tags=None):
    '''
    This function aims to be the gateway for setting up a set of matched-filter
    jobs in a workflow. This function is intended to support multiple
    different ways/codes that could be used for doing this. For now the only
    supported sub-module is one that runs the matched-filtering by setting up
    a serious of matched-filtering jobs, from one executable, to create
    matched-filter triggers covering the full range of science times for which
    there is data and a template bank file.
    Parameters
    -----------
    Workflow : pycbc.workflow.core.Workflow
        The workflow instance that the coincidence jobs will be added to.
    science_segs : ifo-keyed dictionary of ligo.segments.segmentlist instances
        The list of times that are being analysed in this workflow.
    datafind_outs : pycbc.workflow.core.FileList
        An FileList of the datafind files that are needed to obtain the
        data used in the analysis.
    tmplt_banks : pycbc.workflow.core.FileList
        An FileList of the template bank files that will serve as input
        in this stage.
    output_dir : path
        The directory in which output will be stored.
    injection_file : pycbc.workflow.core.File, optional (default=None)
        If given the file containing the simulation file to be sent to these
        jobs on the command line. If not given no file will be sent.
    tags : list of strings (optional, default = [])
        A list of the tagging strings that will be used for all jobs created
        by this call to the workflow. An example might be ['BNSINJECTIONS'] or
        ['NOINJECTIONANALYSIS']. This will be used in output names.
    Returns
    -------
    inspiral_outs : pycbc.workflow.core.FileList
        A list of output files written by this stage. This *will not* contain
        any intermediate products produced within this stage of the workflow.
        If you require access to any intermediate products produced at this
        stage you can call the various sub-functions directly.
    '''
    if tags is None:
        tags = []
    logging.info("Entering matched-filtering setup module.")
    make_analysis_dir(output_dir)
    cp = workflow.cp

    # Parse for options in .ini file
    mfltrMethod = cp.get_opt_tags("workflow-matchedfilter", "matchedfilter-method",
                                  tags)
    

    # Could have a number of choices here
    if mfltrMethod == "WORKFLOW_INDEPENDENT_IFOS":
        logging.info("Adding matched-filter jobs to workflow.")
        if cp.has_option_tags("workflow-matchedfilter",
                              "matchedfilter-link-to-tmpltbank", tags):
            if not cp.has_option_tags("workflow-tmpltbank",
                              "tmpltbank-link-to-matchedfilter", tags):
                errMsg = "If using matchedfilter-link-to-tmpltbank, you should "
                errMsg += "also use tmpltbank-link-to-matchedfilter."
                logging.warn(errMsg)
            linkToTmpltbank = True
        else:
            linkToTmpltbank = False
        if cp.has_option_tags("workflow-matchedfilter",
                              "matchedfilter-compatibility-mode", tags):
            if not linkToTmpltbank:
                errMsg = "Compatibility mode requires that the "
                errMsg += "matchedfilter-link-to-tmpltbank option is also set."
                raise ValueError(errMsg)
            if not cp.has_option_tags("workflow-tmpltbank",
                              "tmpltbank-compatibility-mode", tags):
                errMsg = "If using compatibility mode it must be set both in "
                errMsg += "the template bank and matched-filtering stages."
                raise ValueError(errMsg)
            compatibility_mode = True
        else:
            compatibility_mode = False
        
        # needs to make changes here--> kanchan
        inspiral_outs = setup_matchedfltr_dax_generated(workflow, science_segs,
                                      datafind_outs, tmplt_banks, output_dir,
                                      injection_file=injection_file,
                                      tags=tags,
                                      link_to_tmpltbank=linkToTmpltbank,
                                      compatibility_mode=compatibility_mode)
        
    elif mfltrMethod == "WORKFLOW_MULTIPLE_IFOS":
        logging.info("Adding matched-filter jobs to workflow.")
        inspiral_outs = setup_matchedfltr_dax_generated_multi(workflow,
                                      science_segs, datafind_outs, tmplt_banks,
                                      output_dir, injection_file=injection_file,
                                      tags=tags)
    else:
        errMsg = "Matched filter method not recognized. Must be one of "
        errMsg += "WORKFLOW_INDEPENDENT_IFOS (currently only one option)."
        raise ValueError(errMsg)

    logging.info("Leaving matched-filtering setup module.")
    return inspiral_outs
