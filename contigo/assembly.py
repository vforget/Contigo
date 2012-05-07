#!/usr/bin/python
# Filename: assembly.py

import sys
import re
import math 
# Classes to store assembly data.  Parsers should be written to populate these classes.

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

# Read name nomenclature

read_re = { 'roche': re.compile("^(?P<base_name>[A-Z0-9]+)(?P<paired_end>_(left|right))?(?P<ali_seg>\.[0-9]+\-[0-9]+)?(?P<other_contig>\.(to|fm)[0-9]+)?(?P<paired_contig>\.pr[0-9]+)?"),
            'illumina': re.compile("^(?P<base_name>.*)"),
            'generic': re.compile("^(?P<base_name>.*)")
            }

class RocheReadNameTokens:
    def __init__(self, name, platform):
        self.name = name
        self.base_name = None # String
        self.is_paired = False # Bool
        self.paired_end = ''
        self.aligned_segment = None # Tuple
        self.other_contig = None # String
        self.paired_contig = None # String
        self.tokenize_roche_read_name(platform)
        
    
    def tokenize_roche_read_name(self, platform):
        m = read_re[platform].match(self.name)
        if m:
            self.base_name = m.group('base_name')
            if platform == 'roche':
                self.is_paired = (m.group('paired_end') != None)
                if self.is_paired:
                    self.paired_end = m.group('paired_end')
                if m.group('ali_seg'): 
                    self.aligned_segment = tuple(m.group('ali_seg')[1:].split('-'))
                if m.group('other_contig'):
                    self.other_contig = 'contig' + \
                                        ("%05d" % (int(m.group('other_contig')[3:]),))
                if m.group('paired_contig'):
                    self.paired_contig = 'contig' + \
                                         ("%05d" % (int(m.group('paired_contig')[3:]),))


    
class Assembly:
    ''' Simple class to store assembly i.e. contigs, and some basic assembly statistics'''
    def __init__(self):
        self.contig_statistics = []
        self.qual_histogram_values = None
        self.depth_histogram_values = None
        self.read_length_histogram_values = None
        self.paired_depth_histogram_values = None
        self.avg_insert_size = None
        self.stddev_insert_size = None
        self.num_templates = 0
        self.BothUnmapped = 0
        self.FalsePair = 0
        self.Link = 0
        self.MultiplyMapped = 0
        self.OneUnmapped = 0
        self.SameContig = 0
        self.pxbp = None
        
    def square(self, x):
        return x * x
    
    def update_pair_statistics(self, pair_status):
        self.num_templates = len(pair_status)
        sizes = []
        for pair in pair_status.values():
            self.__dict__[pair.status] += 1
            if pair.status == 'SameContig':
                sizes.append(pair.distance)
        try:
            self.avg_insert_size = float(sum(sizes))/len(sizes)
            self.stddev_insert_size = math.sqrt(sum([self.square(s-self.avg_insert_size) \
                                                         for s in sizes]) / len(sizes))
        except ZeroDivisionError: 
            self.avg_insert_size = 0
            self.stddev_insert_size = 0

    def max(self, attr):
        return max([x.__dict__[attr] for x in self.contig_statistics])
    
    def min(self, attr):
        return min([x.__dict__[attr] for x in self.contig_statistics])
        
    def sum(self,attr):
        return sum([x.__dict__[attr] for x in self.contig_statistics])
        
    def num_contigs(self):
        return len(self.contig_statistics)
    
    def avg_depth(self):
        total_bases_sequenced = 0
        for k in self.depth_histogram():
            total_bases_sequenced += k * self.depth_histogram()[k]
        
        return float(total_bases_sequenced)/self.sum("length")

    def depth_range(self):
        y = [x.avg_depth() for x in self.contig_statistics]
        return min(y), max(y)

    def avg_paired_depth(self):
        total_bases_sequenced = 0
        for k in self.paired_depth_histogram():
            total_bases_sequenced += k * self.paired_depth_histogram()[k]
        
        return float(total_bases_sequenced)/self.sum("length")

    def avg_ins_size(self):
        s = 0.0
        n = 0.0
        for cs in self.contig_statistics:
            s += sum(cs.ins_sizes)
            n += len(cs.ins_sizes)
        try:
            return s/n
        except ZeroDivisionError:
            return 0

    def num_paired(self):
        s = 0
        for cs in self.contig_statistics:
            s += cs.paired_internal
        assert (s % 2) == 0, "Inconsistent number of paired ends."
        return s/2
        
    def paired_depth_range(self):
        y = [x.avg_paired_depth() for x in self.contig_statistics]
        return min(y), max(y)

    def avg_gc(self):
        n, gc = 0, 0
        for cs in self.contig_statistics:
            gc += sum([cs.base_counts["G"], cs.base_counts["C"]])
            n += sum(cs.base_counts.values())
        return (float(gc)/n)*100
    
    def avg_read_length(self):
        c, s, v = 0, 0, None
        for k in self.read_length_histogram():
            v = self.read_length_histogram()[k]
            c += v
            s += v * k
        return float(s)/c

    def median_read_length(self):
        num_reads = self.sum("num_reads")
        sum_of_reads = 0
        for k in self.read_length_histogram():
            sum_of_reads += self.read_length_histogram()[k]
            if sum_of_reads > (num_reads/2):
                return k
        
    def perc_lowqual_range(self):
        y = [x.perc_lowqual() for x in self.contig_statistics]
        return min(y), max(y)

    
    def gc_range(self):
        y = [x.perc_gc() for x in self.contig_statistics]
        return min(y), max(y)
    
    def qual_range(self):
        y = []
        for x in self.contig_statistics: y += x.qual_histogram.keys() 
        return min(y), max(y)
    
    def read_length_range(self):
        ks = self.read_length_histogram().keys()
        return min(ks), max(ks)
    
    def qual_histogram(self):
        if not self.qual_histogram_values:
            self.qual_histogram_values = {}
            for cs in self.contig_statistics:
                for k in cs.qual_histogram:
                    self.qual_histogram_values[k] = \
                        self.qual_histogram_values.get(k,0) + \
                        cs.qual_histogram[k]
        return self.qual_histogram_values
    
    def depth_histogram(self):
        if not self.depth_histogram_values:
            self.depth_histogram_values = {}
            for cs in self.contig_statistics:
                for k in cs.depth_histogram:
                    self.depth_histogram_values[k] = \
                        self.depth_histogram_values.get(k,0) + \
                        cs.depth_histogram[k]
        return self.depth_histogram_values
    
    def read_length_histogram(self):
        if not self.read_length_histogram_values:
            self.read_length_histogram_values = {};
            for cs in self.contig_statistics:
                for k in cs.read_length_histogram:
                    self.read_length_histogram_values[k] = \
                        self.read_length_histogram_values.get(k,0) + \
                        cs.read_length_histogram[k]
        
        return self.read_length_histogram_values

    def paired_depth_histogram(self):
        if not self.paired_depth_histogram_values:
            self.paired_depth_histogram_values = {}
            for cs in self.contig_statistics:
                for k in cs.paired_depth_histogram:
                    self.paired_depth_histogram_values[k] = \
                        self.paired_depth_histogram_values.get(k,0) + \
                        cs.paired_depth_histogram[k]
        return self.paired_depth_histogram_values
    
    
    
    
    
class ContigStatistic:
    
    def __init__(self, name, length, num_reads, base_counts, \
                     depth_histogram, qual_histogram, low_qual_count, \
                     read_length_histogram):
        self.name = name
        self.length = length
        self.num_reads = num_reads
        self.base_counts = base_counts
        self.depth_histogram = depth_histogram
        self.qual_histogram = qual_histogram
        self.low_qual_count = low_qual_count
        self.read_length_histogram = read_length_histogram
        self.wgs = None
        self.paired_internal = None
        self.paired_gap = None,
        self.paired_single = None
        self.ins_sizes = None
        self.paired_depth = None
        self.avg_ins_size_val = 0
        self.avg_paired_depth_val = 0
        self.num_segments = None
        
    def perc_gc(self):
        return (float(self.base_counts["G"] + self.base_counts["C"]) / \
                    sum(self.base_counts.values()))*100
    
    def perc_lowqual(self):
        return (float(self.low_qual_count)/self.length)*100
    
    def avg_depth(self):
        tot = 0
        for k in self.depth_histogram: tot += k * self.depth_histogram[k]
        return float(tot)/self.length

    def avg_ins_size(self):
        try:
            self.avg_ins_size_val = float(sum(self.ins_sizes))/len(self.ins_sizes)
            return self.avg_ins_size_val
        except ZeroDivisionError:
            return self.avg_ins_size_val
        
    def avg_paired_depth(self):
        self.avg_paired_depth_val = float(sum(self.paired_depth))/len(self.paired_depth)
        return self.avg_paired_depth_val

    def med_ins_size(self):
        try:
            self.ins_sizes.sort()
            return self.ins_sizes[len(self.ins_sizes)/2]
        except IndexError:
            return 0
        
    def med_paired_depth(self):
        self.paired_depth.sort()
        return self.paired_depth[len(self.paired_depth)/2]

    def add_read_stats(self, wgs, paired_internal, paired_gap, paired_single, ins_sizes, paired_depth, paired_depth_histogram, templates):
        self.wgs = wgs
        self.paired_internal = paired_internal
        self.paired_gap = paired_gap
        self.paired_single = paired_single
        self.ins_sizes = ins_sizes
        self.paired_depth = paired_depth
        self.paired_depth_histogram = paired_depth_histogram
        self.templates = templates
class Contig:
    ''' Stores contig data i.e.  sequence, qualities, array of reads.'''
    
    BEST_QUAL = { "roche": 64, "illumina": 64, "generic": 64 }
    def __init__(self, name, bases, num_reads, orientation):
        # default start at base 0, but possibly > 0 if reads overhang contig edge
        self.start = 0 
        self.name = name
        self.bases = int(bases)
        self.end = self.start + self.bases
        self.num_reads = int(num_reads)
        self.orientation = orientation
        self.consensus = ""
        self.reads = {}
        self.quality = []
        # Values populated by methods below
        self.depth_values = None
        self.depth_histogram_values = None
        self.base_count_values = None
        self.fastq_value = None
        self.low_qual_count_value = None
        self.qual_histogram_values = None
        self.read_length_histogram_values = None
        
        
    def length(self):
        return len(self.quality)
        
    def to_fastq(self, q):
        '''Convert quality integer to ASCII char.'''
        if q > 93: q += 93
        return chr(q + 33)

    def unpadded_consensus(self):
        ''' Remove pads from consensus '''
        return self.consensus.replace('*', '')
    
    def depth(self):
        '''Calculate depth from read positions.'''
        if not self.depth_values:
            contig_length = len(self.consensus)
            depth_values = [0] * contig_length
            for rn in self.reads:
                rn_start = self.reads[rn].padded_start - 1
                rn_end = self.reads[rn].padded_start + \
                    self.reads[rn].num_padded_bases - 1
                if rn_start < 1: rn_start = 0
                if rn_end > contig_length: rn_end = contig_length
                for i in range(rn_start, rn_end):
                    depth_values[i] += 1
            self.depth_values = []
            for i in range(0, contig_length):
                if self.consensus[i] != '*':
                    self.depth_values.append(depth_values[i])
        return self.depth_values
    
    def read_length_histogram(self):
        if not self.read_length_histogram_values:
            self.read_length_histogram_values = {}
            for r in self.reads:
                rl = self.reads[r].length()
                self.read_length_histogram_values[rl] = \
                    self.read_length_histogram_values.get(rl,0) + 1
        return self.read_length_histogram_values
    

    def depth_histogram(self):
        if not self.depth_histogram_values:
            self.depth_histogram_values = {}
            for x in self.depth():
                self.depth_histogram_values[x] = \
                    self.depth_histogram_values.get(x,0) + 1
    
        return self.depth_histogram_values

    def update_quality_attr(self):
        self.fastq_value = ''
        self.qual_histogram_values = {}
        self.low_qual_count_value = 0
        for x in self.quality:
            self.qual_histogram_values[x] = \
                self.qual_histogram_values.get(x,0) + 1
            self.fastq_value += self.to_fastq(x)
            if x < self.BEST_QUAL["roche"]: self.low_qual_count_value += 1

    def fastq(self):
        if not self.fastq_value:
            self.update_quality_attr()
        return self.fastq_value

    def low_qual_count(self):
        if not self.low_qual_count_value:
            self.update_quality_attr()
        return self.low_qual_count_value
    
    def qual_histogram(self):
        if not self.qual_histogram_values:
            self.update_quality_attr()
        return self.qual_histogram_values
    
    def read_stats(self):
        base_names = set([read.base_name for read in self.reads.values()])
        wgs = 0
        paired_internal = 0
        paired_single = 0
        paired_gap = 0
        ins_sizes = []
        paired_depth = [0] * len(self.consensus)
        templates = {}
        for base_name in base_names:
            ln = base_name + "_left"
            if self.reads.get(ln):
                rn = base_name + "_right"
                if self.reads.get(rn):
                    if (self.reads[ln].paired_contig == self.name and \
                       self.reads[ln].paired_contig == self.reads[rn].paired_contig):
                        paired_internal += 2
                        templates[base_name] = "paired_internal"
                        x = [self.reads[ln].padded_start, self.reads[ln].padded_end,
                             self.reads[rn].padded_start, self.reads[rn].padded_end,
                             ]
                        s, e = min(x), max(x)
                        size = e - s
                        for i in range(s,e):
                            try:
                                paired_depth[i] += 1
                            except IndexError:
                                pass
                        ins_sizes.append(size)
                        
                    elif (self.reads[ln].paired_contig == None and \
                        self.reads[ln].paired_contig == self.reads[rn].paired_contig):
                        templates[base_name] = "paired_internal"
                        paired_internal += 2
                    elif self.reads[ln].paired_contig != self.reads[rn].paired_contig:
                        templates[base_name] = "paired_gap"
                        paired_gap += 2
                    else:
                        print self.reads[ln].name, self.reads[rn].name
                else:
                    templates[base_name] = "paired_single"
                    paired_single += 1
            elif self.reads.get(base_name + "_right"):
                templates[base_name] = "paired_single"
                paired_single += 1
            else:
                wgs += 1
                
        paired_depth_histogram = {}
        for x in paired_depth:
            paired_depth_histogram[x] = \
                paired_depth_histogram.get(x,0) + 1
        
        return (wgs, paired_internal, paired_gap, paired_single, ins_sizes, paired_depth, paired_depth_histogram, templates)
                
                
    
            
    def base_counts(self):
        if not self.base_count_values:
            valid = ("A", "C", "G", "T")
            self.base_count_values = {"A": 0, "C": 0, "G": 0, "T": 0}
            for x in self.consensus:
                if x.upper() in valid: self.base_count_values[x.upper()] += 1

        return self.base_count_values
        

class Read:
    ''' Stores read data i.e. sequence, qualities, contig position.'''
    
    def __init__(self, name, orientation, padded_start, platform):
        self.name = name
        self.read_name = False
        self.orientation = orientation
        self.padded_start = int(padded_start)
        self.qual_clip_start = False
        self.qual_clip_end = False
        self.align_clip_start = False
        self.align_clip_end = False
        self.num_padded_bases = False
        self.sequence = False
        self.length_value = None
        self.padded_end = None
        # 454 instance vars ... should sub-class in future
        self.base_name = None # String
        self.is_paired = False # Bool
        self.paired_end = ''
        self.aligned_segment = None # Tuple
        self.other_contig = None # String
        self.paired_contig = None # String
        self.tokenize_roche_read_name(platform)
        
    def length(self):
        if not self.length_value:
            self.length_value = len(self.sequence.replace('*',''))
        return self.length_value

    def padded_length(self):
        return len(self.sequence)
        

    def tokenize_roche_read_name(self, platform):
        m = read_re[platform].match(self.name)
        if m:
            self.base_name = m.group('base_name')
            if platform == 'roche':
                self.is_paired = (m.group('paired_end') != None)
                if self.is_paired:
                    self.paired_end = m.group('paired_end')
                if m.group('ali_seg'): 
                    self.aligned_segment = tuple(m.group('ali_seg')[1:].split('-'))
                if m.group('other_contig'):
                    self.other_contig = 'contig' + \
                                        ("%05d" % (int(m.group('other_contig')[3:]),))
                if m.group('paired_contig'):
                    self.paired_contig = 'contig' + \
                                         ("%05d" % (int(m.group('paired_contig')[3:]),))
            

            
# Classes and functions to parse and store data from an AGP formatted file.

class Component:
    ''' class for component line '''
    def __init__(self, object_id,  object_beg, object_end, part_number, component_type, component_id, component_beg,  component_end, orientation):
        self.object_id = object_id
        self.object_beg = int(object_beg)
        self.object_end = int(object_end)
        self.part_number = int(part_number)
        self.component_type = component_type
        self.component_id = component_id
        self.component_beg = int(component_beg)
        self.component_end = int(component_end)
        self.orientation = orientation

class Gap:
    ''' class for gap line '''
    def __init__(self, object_id,  object_beg, object_end, part_number, component_type, gap_length, gap_type, linkage):
        self.object_id = object_id
        self.object_beg = int(object_beg)
        self.object_end = int(object_end)
        self.part_number = int(part_number)
        self.component_type = component_type
        self.gap_length = int(gap_length)
        self.gap_type = gap_type
        self.linkage = linkage
        
def load_scaffold_file(fh):
    ''' Parses AGP format.
    Returns hash key\'ed on scaffold name. '''

    scaffolds = {}

    for line in fh:
        columns = line.split()
        component_type = columns[4]
        object_id = columns[0]
        if not scaffolds.get(object_id): scaffolds[object_id] = []
        if component_type != "N": # component
            scaffolds[object_id].append(Component(*columns))
        else: # gap
            scaffolds[object_id].append(Gap(*columns))
    return scaffolds


# Functions to load ACE formatted file from 454

def tokenize_readname(read_name):
    ''' Tokenize read name '''
    


def get_seq(f):
    ''' Get sequence lines.
    Returns string of DNA sequence.'''
    
    indata = True
    data = ""
    seq_re = re.compile('[ACGTQWSMKRYBDHVNXacgtqwsmkrybdhvnx\*]+')
    while indata:
        line = f.readline()
        if re.match(seq_re, line):
            data += line.strip()
        else:
            indata = False
    return data.upper()
        
def get_qual(f):
    ''' Get quality values.
    Returns array of quality values.'''
    
    indata = True
    data = ""
    qual_re = re.compile('\s*[0-9]+')
    while indata:
        line = f.readline()
        if re.match(qual_re, line):
            data += line
        else:
            indata = False
    return [int(x) for x in data.split()]
 
def add_qa_data(data, read):
    '''Parse QA line.
    Updates read.'''

    read.qual_clip_start = int(data[1])
    read.qual_clip_end = int(data[2])
    read.align_clip_start = int(data[3])
    read.align_clip_end = int(data[4])


def load_ace(filename="454Contigs.ace", platform='roche'):
    '''Parses ace file.
    Returns hash keyed on contig name.'''

    f = file(filename, "r")
    line = f.readline()
    data = line.split()
    num_contigs = int(data[1])
    num_reads = int(data[2])
    
    ace = Assembly()
    contig = False
    cc = 0
    rc = 0
    
    while True:
        line = f.readline()
        if len(line) == 0:
            break
    
        if re.match("CO ", line):
            cc += 1
            if contig:
                yield contig
                
            co_data = line.split()
            contig = Contig(co_data[1], co_data[2],
                                     co_data[3], co_data[-1])
            contig.consensus = get_seq(f)
            continue

        if re.match("BQ", line):
            contig.quality = get_qual(f)
            continue

        if re.match("AF ", line):
            rc += 1
            af_data = line.split()
            read = Read(af_data[1], af_data[2], af_data[3], platform)
            contig.reads[read.base_name + read.paired_end] = read
            continue
        
        if re.match("RD ", line):
            rd_data = line.split()
            rt = RocheReadNameTokens(rd_data[1], platform)
            read_name = rt.base_name + rt.paired_end
            contig.reads[read_name].num_padded_bases = int(rd_data[2])
            contig.reads[read_name].padded_end = \
                   contig.reads[read_name].padded_start + \
                   contig.reads[read_name].num_padded_bases
            contig.reads[read_name].sequence = get_seq(f)
            continue
        
        if re.match("QA ", line):
            qa_data = line.split()
            add_qa_data(qa_data, contig.reads[read_name])
            continue
    
    yield contig



def make_scaffold_table_row(sid, item):
    
    rd = None
    r = re.compile('(scaffold|contig)[0]*')
    is_scaffold = False
    if item.component_type == "W":
        link = item.component_id # r.sub('ctg', item.component_id)
        
        rd = [sid,
              # r.sub('scf', item.object_id),
              item.object_id,
              item.object_beg,
              item.object_end,
              (item.object_end - item.object_beg),
              link,
              item.orientation.replace('+', 'U').replace('-', 'C')
              ]
        is_scaffold = True
    elif item.component_type == "N":
        rd = [sid,
              # r.sub('scf', item.object_id),
              item.object_id,
              item.object_beg,
              item.object_end,
              item.gap_length,
              "gap",
              'G'
              ]
    if not rd: exit("Component type \"%s\" not supported" % (item.component_type,))
    return rd, is_scaffold

def make_scaffold_table(options):

    scaffolds = load_scaffold_file(\
        open(options.input_dir + "/" + options.scaffold))
    sid = 0
    sorted_names = scaffolds.keys()
    sorted_names.sort()
    rows = []
    num_scaffolds = 0
    for name in sorted_names: 
        scaffold = scaffolds[name]
        for item in scaffold:
            sid += 1
            row, is_scaffold = \
                make_scaffold_table_row(sid, item)
            rows.append(row)
            if is_scaffold: num_scaffolds += 1
            
    return rows, num_scaffolds

if __name__ == '__main__':
    fh = open(sys.argv[1])
    res = load_454PairStatus(fh)
    for p in res:
        print res[p].distance
    
    
    
