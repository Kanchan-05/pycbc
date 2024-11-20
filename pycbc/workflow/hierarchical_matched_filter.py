from __future__ import division

import os, logging
import math
from pycbc.workflow.core import FileList, make_analysis_dir, File, Node
from ligo import segments
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

        if self.curr_seg_length == self.data_length:
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

    def get_data_times_for_job(self, num_job):
        """ Get the data that this job will read in. """
        # small factor of 0.0001 to avoid float round offs causing us to
        # miss a second at end of segments.
        shift_dur = self.curr_seg[0] + int(self.job_time_shift * num_job\
                                           + 0.0001)
        job_data_seg = self.data_chunk.shift(shift_dur)
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
    
def identify_needed_data(curr_exe_job):
    """ This function will identify the length of data that a specific
    executable needs to analyse and what part of that data is valid (ie.
    inspiral doesn't analyse the first or last 8s of data it reads in).

    Parameters
    -----------
    curr_exe_job : Job
        An instance of the Job class that has a get_valid times method.

    Returns
    --------
    dataLength : float
        The amount of data (in seconds) that each instance of the job must read
        in.
    valid_chunk : ligo.segments.segment
        The times within dataLength for which that jobs output **can** be
        valid (ie. for inspiral this is (72, dataLength-72) as, for a standard
        setup the inspiral job cannot look for triggers in the first 72 or
        last 72 seconds of data read in.)
    valid_length : float
        The maximum length of data each job can be valid for. This is
        abs(valid_segment).
    """
    # Set up the condorJob class for the current executable
    data_lengths, valid_chunks = curr_exe_job.get_valid_times()

    # Begin by getting analysis start and end, and start and end of time
    # that the output file is valid for
    valid_lengths = [abs(valid_chunk) for valid_chunk in valid_chunks]

    return data_lengths, valid_chunks, valid_lengths

# Copied from jobsetup.py
def sngl_ifo_job_setup(workflow, ifo, out_files, curr_exe_job, science_segs,
                       datafind_outs, parents=None,
                       allow_overlap=True):
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
    allow_overlap : boolean (optional, kwarg, default = True)
        If this is set the times that jobs are valid for will be allowed to
        overlap. This may be desired for template banks which may have some
        overlap in the times they cover. This may not be desired for inspiral
        jobs, where you probably want triggers recorded by jobs to not overlap
        at all.

    Returns
    --------
    out_files : pycbc.workflow.core.FileList
        A list of the files that will be generated by this step in the
        workflow.
    """

    ########### (1) ############
    # Get the times that can be analysed and needed data lengths
    data_length, valid_chunk, valid_length = identify_needed_data(curr_exe_job)
    
    exe_tags = curr_exe_job.tags
    # Loop over science segments and set up jobs
    for curr_seg in science_segs:
        ########### (2) ############
        # Initialize the class that identifies how many jobs are needed and the
        # shift between them.
        segmenter = JobSegmenter(data_length, valid_chunk, valid_length,
                                 curr_seg, curr_exe_job)

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

            for parent in curr_parent:
                if len(curr_parent) != 1:
                    bank_tag = [t for t in parent.tags if 'bank' in t.lower()]
                    curr_exe_job.update_current_tags(bank_tag + exe_tags)
                # We should generate unique names automatically, but it is a
                # pain until we can set the output names for all Executables
                node = curr_exe_job.create_node(job_data_seg, job_valid_seg,
                                                parent=parent,
                                                df_parents=curr_dfouts)
                workflow.add_node(node)
                curr_out_files = node.output_files
                # FIXME: Here we remove PSD files if they are coming through.
                #        This should be done in a better way. On to-do list.
                curr_out_files = [i for i in curr_out_files if 'PSD_FILE'\
                                                                 not in i.tags]
                out_files += curr_out_files

    return out_files
# no change till here..

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
#         'pycbc_inspiral'          : PyCBCInspiralExecutable,
#         'pycbc_inspiral_skymax'   : PyCBCInspiralExecutable,
#         'pycbc_multi_inspiral'    : PyCBCMultiInspiralExecutable,
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
                                    tags=None):
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
    # setting matched filtering jobs depending on the tags for hierarchical search 
    logging.info("Setting matched filtering jobs for %s." %(tags))
    
    match_fltr_exe = os.path.basename(cp.get('executables','inspiral'))
    # Select the appropriate class
    exe_class = select_matchedfilter_class(match_fltr_exe)

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
        
        # Modifying sngls job setup function for hierarchical search
        sngl_ifo_job_setup(workflow, ifo, inspiral_outs, job_instance,
                           science_segs[ifo], datafind_outs,
                           parents=tmplt_banks, allow_overlap=False)
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
        # needs to make changes here--> kanchan
        inspiral_outs = setup_matchedfltr_dax_generated(workflow, science_segs,
                                      datafind_outs, tmplt_banks, output_dir,
                                      injection_file=injection_file,
                                      tags=tags)
        
    elif mfltrMethod == "WORKFLOW_MULTIPLE_IFOS":
        logging.info("Adding matched-filter jobs to workflow.")
        inspiral_outs = setup_matchedfltr_dax_generated_multi(workflow,
                                      science_segs, datafind_outs, tmplt_banks,
                                      output_dir, injection_file=injection_file,
                                      tags=tags)
    else:
        errMsg = "Matched filter method not recognized. Must be one of "
        errMsg += "WORKFLOW_INDEPENDENT_IFOS or WORKFLOW_MULTIPLE_IFOS."
        raise ValueError(errMsg)

    logging.info("Leaving matched-filtering setup module.")
    return inspiral_outs