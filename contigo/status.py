#!/usr/bin/python

import os
import sys
import re

READSTATUS_TO_RGB = {
    'Assembled': ((0,0,0), (100,100,100)),
    'Repeat': ((160,80,0), (235,172,0)),
    'PartiallyAssembled': ((165,165,165), (0,0,0)),
    'Outlier': ((165,165,165), (0,0,0)),
    'Singleton': ((165,165,165), (0,0,0)),
    'TooShort': ((165,165,165), (0,0,0)),
    'paired_internal': ((0,0,130), (100,100,255)),
    'paired_gap': ((0,130,0), (100,255,100)),
    'paired_single': ((130,0,0), (255,100,100)) }

READPAIRSTATUS_TO_RGB = {
    'SameContig': ((100,100,255),(0,0,130)),
    'Link': ((0,130,0), (100,255,100)),
    'FalsePair': ((130,0,0), (255,100,100)),
    'BothUnmapped': ((130,0,130), (255,100,255)), 
    'MultiplyMapped': ((130,130,0), (255,255,100)),
    'OneUnmapped': ((165,165,165), (0,0,0)) }


class PairStatus:

    def __init__(self, template, status, distance, left_contig, left_pos, left_dir,
                 right_contig, right_pos, right_dir, left_distance,
                 right_distance):
        self.template = template
        self.status = status
        self.distance = self.sanitize_value(distance)
        self.left_contig = self.sanitize_value(left_contig)
        self.left_pos = self.sanitize_value(left_pos)
        self.left_dir = self.sanitize_value(left_dir)
        self.right_contig = self.sanitize_value(right_contig)
        self.right_pos = self.sanitize_value(right_pos)
        self.right_dir = self.sanitize_value(right_dir)
        self.left_distance = self.sanitize_value(left_distance)
        self.right_distance = self.sanitize_value(right_distance)
        
    def sanitize_value(self, val):
        if val in ('-', '', '\n'):
            return None
        else:
            if re.match("^[0-9]+$", val):
                return int(val)
            else:
                return val

def load_pair_status(fh):
    res = {}
    line = fh.readline()
    fields = line.split('\t')
    if fields[0] != 'Template':
        p = PairStatus(*fields)
        res[p.template] = p
        
    for line in fh:
        fields = line.split('\t')
        p = PairStatus(*fields)
        res[p.template] = p
    return res

def process_pairstatus(options):
    pair_status = {}
    sys.stderr.write("PairStatus ... ")
    if options.pair_status and \
            os.path.exists(options.input_dir + "/" + options.pair_status):
        filename = options.input_dir + "/" + options.pair_status
        pair_status = load_pair_status(open(filename))
    else:
        sys.stderr.write("FILE: %s does not exist ... skipping ... " % options.pair_status)
    sys.stderr.write("done.\n")
    return pair_status

def load_read_status(filename="454ReadStatus.txt"):
    ''' Parses read status file.
    Returns hash key\'ed on read name.'''

    head_re = re.compile('Accno\tRead\sStatus')
    f = file(filename, "r")
    read_status = {}
    line = f.readline()
    if not re.match(head_re,line):
        name, status = re.match("(\S+)\s+(\S+)", line).groups()
        try:
            read_status[name] = READSTATUS_TO_RGB[status];
        except KeyError:
            sys.stderr.write(
                "No color code for %s (%s), defaulting to Assembled\n" % (status, name))
            read_status[name] = READSTATUS_TO_RGB['Assembled'];

    for line in f:
        name, status = re.match("(\S+)\s+(\S+)", line).groups()
        try:
            read_status[name] = READSTATUS_TO_RGB[status];
        except KeyError:
            sys.stderr.write(
                "No color code for %s (%s), defaulting to Assembled\n" % (status, name))
            read_status[name] = READSTATUS_TO_RGB['Assembled'];
    f.close()
    return read_status


def process_readstatus(options):
    read_status = {}
    sys.stderr.write("ReadStatus ... ")
    if options.read_status and \
            os.path.exists(options.input_dir + "/" + options.read_status):
        filename = options.input_dir + "/" + options.read_status
        read_status = load_read_status(filename)
    else:
        sys.stderr.write("FILE: %s does not exist ... skipping ... " % options.read_status)
    sys.stderr.write("done.\n")
    return read_status

