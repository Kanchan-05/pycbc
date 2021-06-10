# Copyright (C) 2013  Ian Harry
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

#import logging
#import math, os
import os
import lal
from ligo import segments
#import Pegasus.DAX3 as dax
from pycbc.workflow.core import File, FileList, Node
from hierarchical_core import Executable


def int_gps_time_to_str(t):
    """Takes an integer GPS time, either given as int or lal.LIGOTimeGPS, and
    converts it to a string. If a LIGOTimeGPS with nonzero decimal part is
    given, raises a ValueError."""

    if isinstance(t, int):
        return str(t)
    elif isinstance(t, float):
        # Wouldn't this just work generically?
        int_t = int(t)
        if abs(t - int_t) > 0.:
            raise ValueError('Need an integer GPS time, got %s' % str(t))
        return str(int_t)
    elif isinstance(t, lal.LIGOTimeGPS):
        if t.gpsNanoSeconds == 0:
            return str(t.gpsSeconds)
        else:
            raise ValueError('Need an integer GPS time, got %s' % str(t))
    else:
        err_msg = "Didn't understand input type {}".format(type(t))
        raise ValueError(err_msg)
        
# For future: add PyCBCHierarchicalInspiralExecutable to the jobsetup only
class PyCBCHierarchicalInspiralExecutable(Executable):
    """ The class used to create jobs for pycbc_hierarchical_inspiral Executable. """

    current_retention_level = Executable.ALL_TRIGGERS
    file_input_options = ['--gating-file','--connection-file','--coarse-bank-file','--fine-bank-file','--coarse-trigger-files']
    
    def __init__(self, cp, exe_name, ifo=None, out_dir=None,
                 injection_file=None, tags=None, reuse_executable=False):
        
        if tags is None:
            tags = []
        # compatibilities with v1.16.13    
        super(PyCBCHierarchicalInspiralExecutable, self).__init__(cp, exe_name, None, ifo,out_dir, tags=tags, reuse_executable=reuse_executable)
        self.cp = cp
        self.set_memory(2000)
        self.injection_file = injection_file
        self.ext = 'stage2.hdf'

        self.num_threads = 1
        if self.get_opt('processing-scheme') is not None:
            stxt = self.get_opt('processing-scheme')
            if len(stxt.split(':')) > 1:
                self.num_threads = stxt.split(':')[1]
          
    # modifications specific to hierarchical 
    def create_node(self, data_seg, valid_seg, parent=None, dfParents=None, tags=None):
        if tags is None:
            tags = []
        node = Node(self, valid_seg=valid_seg)
        if not self.has_opt('pad-data'):
            raise ValueError("The option pad-data is a required option of "
                             "%s. Please check the ini file." % self.name)
        pad_data = int(self.get_opt('pad-data'))

        # set remaining options flags
        node.add_opt('--gps-start-time',
                     int_gps_time_to_str(data_seg[0] + pad_data))
        node.add_opt('--gps-end-time',
                     int_gps_time_to_str(data_seg[1] - pad_data))
        node.add_opt('--trig-start-time', int_gps_time_to_str(valid_seg[0]))
        node.add_opt('--trig-end-time', int_gps_time_to_str(valid_seg[1]))
        node.add_profile('condor', 'request_cpus', self.num_threads)
        
        if self.injection_file is not None:
            node.add_input_opt('--injection-file', self.injection_file)
        
        # set the input and output files
        fil = node.new_output_file_opt(valid_seg, self.ext, '--output', tags=tags,
                         store_file=self.retain_files, use_tmp_subdirs=True)
        fil.add_metadata('data_seg', data_seg)
        node.add_input_opt('--bank-file', parent)
        if dfParents is not None:
            node.add_input_list_opt('--frame-files', dfParents)

        return node

    def get_valid_times(self):
        """ Determine possible dimensions of needed input and valid output
        """

        if self.cp.has_option('workflow-matchedfilter',
                                                      'min-analysis-segments'):
            min_analysis_segs = int(self.cp.get('workflow-matchedfilter',
                                                'min-analysis-segments'))
        else:
            min_analysis_segs = 0

        if self.cp.has_option('workflow-matchedfilter',
                                                      'max-analysis-segments'):
            max_analysis_segs = int(self.cp.get('workflow-matchedfilter',
                                                'max-analysis-segments'))
        else:
            # Choose ridiculously large default value
            max_analysis_segs = 1000

        if self.cp.has_option('workflow-matchedfilter', 'min-analysis-length'):
            min_analysis_length = int(self.cp.get('workflow-matchedfilter',
                                                  'min-analysis-length'))
        else:
            min_analysis_length = 0

        if self.cp.has_option('workflow-matchedfilter', 'max-analysis-length'):
            max_analysis_length = int(self.cp.get('workflow-matchedfilter',
                                                  'max-analysis-length'))
        else:
            # Choose a ridiculously large default value
            max_analysis_length = 100000

        segment_length = int(self.get_opt('segment-length'))
        pad_data = 0
        if self.has_opt('pad-data'):
            pad_data += int(self.get_opt('pad-data'))

        # NOTE: Currently the tapered data is ignored as it is short and
        #       will lie within the segment start/end pad. This means that
        #       the tapered data *will* be used for PSD estimation (but this
        #       effect should be small). It will also be in the data segments
        #       used for SNR generation (when in the middle of a data segment
        #       where zero-padding is not being used) but the templates should
        #       not be long enough to use this data assuming segment start/end
        #       pad take normal values. When using zero-padding this data will
        #       be used for SNR generation.

        #if self.has_opt('taper-data'):
        #    pad_data += int(self.get_opt( 'taper-data' ))
        if self.has_opt('allow-zero-padding'):
            self.zero_padding=True
        else:
            self.zero_padding=False

        start_pad = int(self.get_opt( 'segment-start-pad'))
        end_pad = int(self.get_opt('segment-end-pad'))

        seg_ranges = range(min_analysis_segs, max_analysis_segs + 1)
        data_lengths = []
        valid_regions = []
        for nsegs in seg_ranges:
            analysis_length = (segment_length - start_pad - end_pad) * nsegs
            if not self.zero_padding:
                data_length = analysis_length + pad_data * 2 \
                              + start_pad + end_pad
                start = pad_data + start_pad
                end = data_length - pad_data - end_pad
            else:
                data_length = analysis_length + pad_data * 2
                start = pad_data
                end = data_length - pad_data
            if data_length > max_analysis_length: continue
            if data_length < min_analysis_length: continue
            data_lengths += [data_length]
            valid_regions += [segments.segment(start, end)]
        # If min_analysis_length is given, ensure that it is added as an option
        # for job analysis length.
        if min_analysis_length:
            data_length = min_analysis_length
            if not self.zero_padding:
                start = pad_data + start_pad
                end = data_length - pad_data - end_pad
            else:
                start = pad_data
                end = data_length - pad_data
            if end > start:
                data_lengths += [data_length]
                valid_regions += [segments.segment(start, end)]

        return data_lengths, valid_regions

    def zero_pad_data_extend(self, job_data_seg, curr_seg):
        """When using zero padding, *all* data is analysable, but the setup
        functions must include the padding data where it is available so that
        we are not zero-padding in the middle of science segments. This
        function takes a job_data_seg, that is chosen for a particular node
        and extends it with segment-start-pad and segment-end-pad if that
        data is available.
        """
        if self.zero_padding is False:
            return job_data_seg
        else:
            start_pad = int(self.get_opt( 'segment-start-pad'))
            end_pad = int(self.get_opt('segment-end-pad'))
            new_data_start = max(curr_seg[0], job_data_seg[0] - start_pad)
            new_data_end = min(curr_seg[1], job_data_seg[1] + end_pad)
            new_data_seg = segments.segment([new_data_start, new_data_end])
            return new_data_seg
        
