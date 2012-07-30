#!/usr/bin/python
# Filename: render.py

import datetime
import sys
import locale
import os

def set_locale():
    version_info = sys.version_info
    print version_info
    if version_info[0] != 2:
        exit("Contigo requires Python version 2.x.")

    if version_info[1] < 7:
        locale.setlocale(locale.LC_ALL, 'en_US')
    if version_info[1] >= 7:
        locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')


HTML_DIR = os.path.join(os.path.dirname(__file__)) + "/static/html/"

def html_header():
    return """ TEXT """

def html_body():
    return """
Manual
Why the name?
Why Yet-Another-Assembly Viewer?
How's it made?
"""


def html_footer():
    return " TEXT "

def interface(of, options, contig_table, scaffold_table):
    html = open(HTML_DIR + "template.html").read()
    html = html.replace("CONTIG_TABLE_HTML", contig_table)
    html = html.replace("SCAFFOLD_TABLE_HTML", scaffold_table)
    html = html.replace("DATE_CREATED", datetime.datetime.now().strftime("%d %b, %Y"))
    html = html.replace("ASSEMBLY_NAME", options.assembly_name)
    print >> of, html
    
def contig_table_header():
    s = '<table id="contig_table">'
    s += """<colgroup>
<col>
<col>
<col>
<col>
<col>
<col>
<col>
<col>
</colgroup>
<tr class="header">
<td>Contig</td>
<td>Length</td>
<td>Read Count</td>
<td>Read Depth</td>
<td>%&lt;Q64</td>
<td>%GC</td>
<!--
<td>Templ.<br/>Depth</td>
<td>Avg.<br/>Insert</td>
-->
        </tr>"""
    return s
    
def contig_table_row(contig_statistic, link):
    set_locale()
    s = "<tr>"
    s += "<td>"
    s += '</td><td>'.join([
        link,
        '<a href="contigo.html" onclick="loadFasta(\'%s\'); return false;" title="View Contig Sequence in FASTA format">%s</a>' % (contig_statistic.name, locale.format("%d", contig_statistic.length, grouping=True),),
        # '<a href="contigo.html" onclick="loadFasta(\'%s\'); return false;" title="View Contig Sequence in FASTA format">%s</a>' % (contig_statistic.name, contig_statistic.length),
        '<a href="contigo.html" onclick="loadReads(\'%s\', %s); return false;" title="View Contig Read Sequences in FASTA format">%s</a>' % (contig_statistic.name, contig_statistic.num_segments, locale.format("%d", contig_statistic.num_reads, grouping=True),),
        ("%0.2f" % contig_statistic.avg_depth()),
        '<a href="contigo.html" onclick="loadQuality(\'%s\'); return false;" title="View Contig Qualities FASTA format">%0.2f</a>' % (contig_statistic.name, contig_statistic.perc_lowqual(),),
        ("%0.2f" % contig_statistic.perc_gc())])
    s += "</td>"
    #s += "<td>%0.2f</td><td>%0.2f</td>" % (contig_statistic.avg_paired_depth(), \
    #                                       contig_statistic.avg_ins_size()/1000
    #                                       )
    s += "</tr>\n"
    return s

def contig_table_footer():
    return "</table>"


def scaffold_table_header():
    s = '<table id="scaffold_table">'
    s += """<colgroup>
<col style="width: 30px;">
<col style="width: 30px;">
<col style="width: 30px;">
<col style="width: 30px;">
</colgroup>
<tr>
<td>Id</td>
<td>Scaffold</td>
<td>Start</td>
<td>End</td>
<td>Len.</td>
<td>Item</td>
<td>Orien.</td>
</tr>"""
    return s

def scaffold_table_row(row_data):
    s = "<tr><td>"
    s += "</td><td>".join([str(x) for x in row_data])
    s += "</td></tr>"
    return s
    
def scaffold_table_footer():
    s = '</table>'
    return s

