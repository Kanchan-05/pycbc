#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
This module is responsible for setting up the segment generation stage of
workflows. For details about this module and its capabilities see here:
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/ahope/segments.html
"""

import sys, os,math,logging
from six.moves import urllib
from six.moves.urllib.request import pathname2url
from six.moves.urllib.parse import urljoin
import numpy
import lal
from ligo import segments
from glue.ligolw import table, lsctables, ligolw
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import segments as ligolw_segments
from glue.ligolw.utils import process as ligolw_process
from pycbc.workflow.core import File
from pycbc.io.ligolw import LIGOLWContentHandler, create_process_table

# class ContentHandler(ligolw.LIGOLWContentHandler):
#     pass

# lsctables.use_in(ContentHandler)

class SegFile(File):
    '''
    This class inherits from the File class, and is designed to store
    workflow output files containing a segment dict. This is identical in
    usage to File except for an additional kwarg for holding the
    segment dictionary, if it is known at workflow run time.
    '''
    def __init__(self, ifo_list, description, valid_segment,
                 segment_dict=None, seg_summ_dict=None, **kwargs):
        """
        See File.__init__ for a full set of documentation for how to
        call this class. The only thing unique and added to this class is
        the optional segment_dict. NOTE that while segment_dict is a
        ligo.segments.segmentlistdict rather than the usual dict[ifo]
        we key by dict[ifo:name].

        Parameters
        ------------
        ifo_list : string or list (required)
            See File.__init__
        description : string (required)
            See File.__init__
        segment : ligo.segments.segment or ligo.segments.segmentlist
            See File.__init__
        segment_dict : ligo.segments.segmentlistdict (optional, default=None)
            A ligo.segments.segmentlistdict covering the times covered by the
            segmentlistdict associated with this file.
            Can be added by setting self.segment_dict after initializing an
            instance of the class.

        """
        super(SegFile, self).__init__(ifo_list, description, valid_segment,
                                      **kwargs)
        # To avoid confusion with the segment_list property of the parent class
        # we refer to this as valid_segments here
        self.valid_segments = self.segment_list
        self.segment_dict = segment_dict
        self.seg_summ_dict = seg_summ_dict

    @classmethod
    def from_segment_list(cls, description, segmentlist, name, ifo,
                          seg_summ_list=None, **kwargs):
        """ Initialize a SegFile object from a segmentlist.

        Parameters
        ------------
        description : string (required)
            See File.__init__
        segmentlist : ligo.segments.segmentslist
            The segment list that will be stored in this file.
        name : str
            The name of the segment lists to be stored in the file.
        ifo : str
            The ifo of the segment lists to be stored in this file.
        seg_summ_list : ligo.segments.segmentslist (OPTIONAL)
            Specify the segment_summary segmentlist that goes along with the
            segmentlist. Default=None, in this case segment_summary is taken
            from the valid_segment of the SegFile class.
        """
        seglistdict = segments.segmentlistdict()
        seglistdict[ifo + ':' + name] = segmentlist
        seg_summ_dict = None
        if seg_summ_list is not None:
            seg_summ_dict = segments.segmentlistdict()
            seg_summ_dict[ifo + ':' + name] = seg_summ_list
        return cls.from_segment_list_dict(description, seglistdict,
                                          seg_summ_dict=seg_summ_dict, **kwargs)

    @classmethod
    def from_multi_segment_list(cls, description, segmentlists, names, ifos,
                                seg_summ_lists=None, **kwargs):
        """ Initialize a SegFile object from a list of segmentlists.

        Parameters
        ------------
        description : string (required)
            See File.__init__
        segmentlists : List of ligo.segments.segmentslist
            List of segment lists that will be stored in this file.
        names : List of str
            List of names of the segment lists to be stored in the file.
        ifos : str
            List of ifos of the segment lists to be stored in this file.
        seg_summ_lists : ligo.segments.segmentslist (OPTIONAL)
            Specify the segment_summary segmentlists that go along with the
            segmentlists. Default=None, in this case segment_summary is taken
            from the valid_segment of the SegFile class.
        """
        seglistdict = segments.segmentlistdict()
        for name, ifo, segmentlist in zip(names, ifos, segmentlists):
            seglistdict[ifo + ':' + name] = segmentlist
        if seg_summ_lists is not None:
            seg_summ_dict = segments.segmentlistdict()
            for name, ifo, seg_summ_list in zip(names, ifos, seg_summ_lists):
                seg_summ_dict[ifo + ':' + name] = seg_summ_list
        else:
            seg_summ_dict = None

        return cls.from_segment_list_dict(description, seglistdict,
                                         seg_summ_dict=seg_summ_dict, **kwargs)

    @classmethod
    def from_segment_list_dict(cls, description, segmentlistdict,
                               ifo_list=None, valid_segment=None,
                               file_exists=False, seg_summ_dict=None,
                               **kwargs):
        """ Initialize a SegFile object from a segmentlistdict.

        Parameters
        ------------
        description : string (required)
            See File.__init__
        segmentlistdict : ligo.segments.segmentslistdict
            See SegFile.__init__
        ifo_list : string or list (optional)
            See File.__init__, if not given a list of all ifos in the
            segmentlistdict object will be used
        valid_segment : ligo.segments.segment or ligo.segments.segmentlist
            See File.__init__, if not given the extent of all segments in the
            segmentlistdict is used.
        file_exists : boolean (default = False)
            If provided and set to True it is assumed that this file already
            exists on disk and so there is no need to write again.
        seg_summ_dict : ligo.segments.segmentslistdict
            Optional. See SegFile.__init__.
        """
        if ifo_list is None:
            ifo_set = set([i.split(':')[0] for i in segmentlistdict.keys()])
            ifo_list = list(ifo_set)
            ifo_list.sort()
        if valid_segment is None:
            if seg_summ_dict and \
                    numpy.any([len(v) for _, v in seg_summ_dict.items()]):
                # Only come here if seg_summ_dict is supplied and it is
                # not empty.
                valid_segment = seg_summ_dict.extent_all()
            else:
                try:
                    valid_segment = segmentlistdict.extent_all()
                except:
                    # Numpty probably didn't supply a
                    # ligo.segments.segmentlistdict
                    segmentlistdict=segments.segmentlistdict(segmentlistdict)
                    try:
                        valid_segment = segmentlistdict.extent_all()
                    except ValueError:
                        # No segment_summary and segment list is empty
                        # Setting valid segment now is hard!
                        warn_msg = "No information with which to set valid "
                        warn_msg += "segment."
                        logging.warn(warn_msg)
                        valid_segment = segments.segment([0,1])
        instnc = cls(ifo_list, description, valid_segment,
                     segment_dict=segmentlistdict, seg_summ_dict=seg_summ_dict,
                     **kwargs)
        if not file_exists:
            instnc.to_segment_xml()
        else:
            instnc.add_pfn(urljoin('file:', pathname2url(instnc.storage_path)),
                           site='local')
        return instnc

    @classmethod
    def from_segment_xml(cls, xml_file, **kwargs):
        """
        Read a ligo.segments.segmentlist from the file object file containing an
        xml segment table.

        Parameters
        -----------
        xml_file : file object
            file object for segment xml file
        """
        # load xmldocument and SegmentDefTable and SegmentTables
        fp = open(xml_file, 'rb')
        xmldoc = ligolw_utils.load_fileobj(
                fp, compress='auto', contenthandler=LIGOLWContentHandler)

        seg_def_table = lsctables.SegmentDefTable.get_table(xmldoc)
        seg_table = lsctables.SegmentTable.get_table(xmldoc)
        seg_sum_table = lsctables.SegmentSumTable.get_table(xmldoc)

        segs = segments.segmentlistdict()
        seg_summ = segments.segmentlistdict()

        seg_id = {}
        for seg_def in seg_def_table:
            # Here we want to encode ifo and segment name
            full_channel_name = ':'.join([str(seg_def.ifos),
                                          str(seg_def.name)])
            seg_id[int(seg_def.segment_def_id)] = full_channel_name
            segs[full_channel_name] = segments.segmentlist()
            seg_summ[full_channel_name] = segments.segmentlist()

        for seg in seg_table:
            seg_obj = segments.segment(
                    lal.LIGOTimeGPS(seg.start_time, seg.start_time_ns),
                    lal.LIGOTimeGPS(seg.end_time, seg.end_time_ns))
            segs[seg_id[int(seg.segment_def_id)]].append(seg_obj)

        for seg in seg_sum_table:
            seg_obj = segments.segment(
                    lal.LIGOTimeGPS(seg.start_time, seg.start_time_ns),
                    lal.LIGOTimeGPS(seg.end_time, seg.end_time_ns))
            seg_summ[seg_id[int(seg.segment_def_id)]].append(seg_obj)

        for seg_name in seg_id.values():
            segs[seg_name] = segs[seg_name].coalesce()

        xmldoc.unlink()
        fp.close()
        curr_url = urllib.parse.urlunparse(['file', 'localhost', xml_file,
                                            None, None, None])

        return cls.from_segment_list_dict('SEGMENTS', segs, file_url=curr_url,
                                          file_exists=True,
                                          seg_summ_dict=seg_summ, **kwargs)

    def remove_short_sci_segs(self, minSegLength):
        """
        Function to remove all science segments
        shorter than a specific length. Also updates the file on disk to remove
        these segments.

        Parameters
        -----------
        minSegLength : int
            Maximum length of science segments. Segments shorter than this will
            be removed.
        """
        newsegment_list = segments.segmentlist()
        for key, seglist in self.segment_dict.items():
            newsegment_list = segments.segmentlist()
            for seg in seglist:
                if abs(seg) > minSegLength:
                    newsegment_list.append(seg)
            newsegment_list.coalesce()
            self.segment_dict[key] = newsegment_list
        self.to_segment_xml(override_file_if_exists=True)

    def return_union_seglist(self):
        return self.segment_dict.union(self.segment_dict.keys())

    def parse_segdict_key(self, key):
        """
        Return ifo and name from the segdict key.
        """
        splt = key.split(':')
        if len(splt) == 2:
            return splt[0], splt[1]
        else:
            err_msg = "Key should be of the format 'ifo:name', got %s." %(key,)
            raise ValueError(err_msg)

    def to_segment_xml(self, override_file_if_exists=False):
        """
        Write the segment list in self.segmentList to self.storage_path.
        """
        # create XML doc and add process table
        outdoc = ligolw.Document()
        outdoc.appendChild(ligolw.LIGO_LW())
        process = create_process_table(outdoc)

        for key, seglist in self.segment_dict.items():
            ifo, name = self.parse_segdict_key(key)
            # Ensure we have LIGOTimeGPS
            fsegs = [(lal.LIGOTimeGPS(seg[0]),
                      lal.LIGOTimeGPS(seg[1])) for seg in seglist]

            if self.seg_summ_dict is None:
                vsegs = [(lal.LIGOTimeGPS(seg[0]),
                          lal.LIGOTimeGPS(seg[1])) \
                         for seg in self.valid_segments]
            else:
                vsegs = [(lal.LIGOTimeGPS(seg[0]),
                          lal.LIGOTimeGPS(seg[1])) \
                         for seg in self.seg_summ_dict[key]]

            # Add using glue library to set all segment tables
            with ligolw_segments.LigolwSegments(outdoc, process) as x:
                x.add(ligolw_segments.LigolwSegmentList(active=fsegs,
                                    instruments=set([ifo]), name=name,
                                    version=1, valid=vsegs))

        # write file
        url = urljoin('file:', pathname2url(self.storage_path))
        if not override_file_if_exists or not self.has_pfn(url, site='local'):
            self.add_pfn(url, site='local')
        ligolw_utils.write_filename(outdoc, self.storage_path)