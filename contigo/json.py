#!/usr/bin/python
# Filename: json.py
# JSON methods for contigo

def assembly_chart_parameters():
    
    return 'CONTIGO.tick_labels = '\
        '{ "length": ["0", "10", "100", "1k", '\
        '"10k", "100k", "1M", "10M"], '\
        '"avg_depth": ["0", "10", "100", "10^3","10^4", "10^5"] };'\
        ' CONTIGO.axis_labels = '\
        '{ "length": "Contig Length", "avg_gc": "GC Content", '\
        '"avg_depth": "Read Depth", "read_length": "Read Length",'\
        '"avg_ins_size": "Ins. Size", "avg_paired_depth": "Tmpl. Depth",'\
        '"perc_low_qual": "Low Quality", "num_reads": "Read Count"'\
        '};'\
        ' CONTIGO.is_log = { "length": true, "avg_depth": true, '\
        '"avg_gc": false, "quality": false, "read_length": false, '\
        '"perc_low_qual": false, "num_reads": false, "avg_ins_size": false, "avg_paired_depth": false, "avg_paired_depth": false'\
        '};'\
        ' CONTIGO.bin_sizes = { "length": 100, "avg_gc": 1, '\
        '"avg_depth": 5, "avg_ins_size": 1, "avg_paired_depth": 1, "perc_low_qual": 1, "num_reads": 5};'


def build_contig_statistic(d):
    s = '{ "name": "%s",' % d.name
    s += '"length": %s,' % d.length
    s += '"avg_depth": %0.0f,' % d.avg_depth()
    s += '"avg_gc": %0.2f,' % d.perc_gc()
    s += '"perc_low_qual": %0.2f,' % d.perc_lowqual()
    s += '"num_reads": %0.2f,' % d.num_reads
    s += '"avg_ins_size": %0.2f,' % (d.avg_ins_size(),)
    s += '"avg_paired_depth": %0.2f,' % (d.avg_paired_depth(),)
    s += '"dzi": "dzi/' + d.name + '.dzi"'
    s += ' }'
    return s
    

def build_assembly_statistic(contigs):
    s = "CONTIGO.assembly = [ "
    s += ",".join([build_contig_statistic(contig) for contig in contigs])
    s += "];"
    return s


def contig_data(contig):
    return "({ 'consensus': '%s',"\
           "'quality': [%s], 'depth': [%s], 'padded_start': %s })\n" % \
           (contig.unpadded_consensus(),  
            ','.join([str(x) for x in contig.quality]),
            ','.join([str(x) for x in contig.depth()]),
            contig.start)


def contig_statistics(asm):
    s = "CONTIGO.num_contigs = %s;" % asm.num_contigs()
    s += "CONTIGO.num_reads = %s;" % asm.sum("num_reads")
    s += "CONTIGO.assembled_bases = %s;" % asm.sum("length")
    s += "CONTIGO.avg_depth = %0.2f;" % asm.avg_depth()
    s += "CONTIGO.avg_paired_depth = %0.1f;" % asm.avg_paired_depth()
    s += "CONTIGO.avg_ins_size = %0.1f;" % (asm.avg_insert_size/1000)
    s += "CONTIGO.stddev_ins_size = %0.1f;" % (asm.stddev_insert_size/1000)
    s += "CONTIGO.num_paired = %s;" % (asm.SameContig)
    s += "CONTIGO.num_templates = %s;" % (asm.num_templates)
    s += "CONTIGO.num_link = %s;" % (asm.Link)
    s += "CONTIGO.num_false = %s;" % (asm.FalsePair)
    s += "CONTIGO.num_unmapped = %s;" % (asm.BothUnmapped)
    s += "CONTIGO.num_multiple = %s;" % (asm.MultiplyMapped)
    s += "CONTIGO.num_half = %s;" % (asm.OneUnmapped)
    s += "CONTIGO.avg_gc = %0.2f;" % asm.avg_gc()
    s += "CONTIGO.avg_read_length = %0.0f;" % asm.avg_read_length()
    s += "CONTIGO.median_read_length = %0.0f;" % asm.median_read_length()
    s += "CONTIGO.ranges = { max: { length: %s, "\
        "avg_depth: %s, avg_gc: %s, quality: %s, read_length: %s, "\
        "avg_ins_size: %s, avg_paired_depth: %s, perc_low_qual: %s, num_reads: %s}, "\
        "min: { length: %s, avg_depth: %s, avg_gc: %s, "\
        "quality: %s, read_length: %s, "\
        "avg_ins_size: %s, avg_paired_depth: %s, perc_low_qual: %s, num_reads: %s"\
        " } };" % (asm.max("length"), 
                   asm.depth_range()[1],
                   asm.gc_range()[1], 
                   asm.qual_range()[1],
                   asm.read_length_range()[1],
                   asm.max("avg_ins_size_val"),
                   asm.max("avg_paired_depth_val"),
                   asm.perc_lowqual_range()[1],
                   asm.max("num_reads"),
                   asm.min("length"), 
                   asm.depth_range()[0], 
                   asm.gc_range()[0], 
                   asm.qual_range()[0], 
                   asm.read_length_range()[0],
                   asm.min("avg_ins_size_val"),
                   asm.min("avg_paired_depth_val"),
                   asm.perc_lowqual_range()[0],
                   asm.min("num_reads")
                   )
                           
    
    s += 'CONTIGO.hist_data = { "avg_depth": [%s], '\
        '"quality": [%s], "read_length": [%s], "avg_paired_depth": [%s] };' % \
        ((",".join([str(asm.depth_histogram().get(i, 0)) \
                        for i in range(0, max(asm.depth_histogram().keys())+1)])),
         (",".join([str(asm.qual_histogram().get(i, 0)) \
                        for i in range(0, max(asm.qual_histogram().keys())+1)])), 
         (",".join([str(asm.read_length_histogram().get(i, 0)) \
                        for i in range(0, max(\
                                 asm.read_length_histogram().keys())+1)])),
         (",".join([str(asm.paired_depth_histogram().get(i, 0)) \
                        for i in range(0, max(asm.paired_depth_histogram().keys())+1)]))

         )
    
         
    s += build_assembly_statistic(asm.contig_statistics)
    s += assembly_chart_parameters()
    return s
