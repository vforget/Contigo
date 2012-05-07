#!/usr/bin/python
# Filename: zoomtig.py
# Copyright Vincenzo Forgetta, 2009-.
# For use with Seadragon AJAX.
# Create tiled segmented image for a contig read assembly.
################################################################################
import Image
import ImageFont
import ImageDraw
# import operator
import sys
import math
import os
import assembly
import multiprocessing



class ImageConfig:
    pass

def calculate(func, args):
    func(*args)
    
def calculatestar(args):
    return calculate(*args)

def get_font_dim(fn):
    return ImageDraw.Draw(
        Image.new("RGB", (1000,1000), (255,255,255))).textsize("A", font=fn)
    
FONT_DIR = sys.path[0] + "/../contigo/static/fonts/"
TTF = FONT_DIR + "Verdana.ttf"


BASE_COLORS = { 
    'A': (255,80,80), 'C': (80,255,80), 'G': (80,80,255), 'T': (255,255,80), 
    'a': (255,80,80), 'c': (80,255,80), 'g': (80,80,255), 't': (255,255,80),
    'n': (255,80,255), 'N': (255,80,255),
    '*': (80,255,255)
    }

SMALL_FONT = ImageFont.load(FONT_DIR + "fixed.6.pil")
# SMALL_FONT = ImageFont.truetype(sys.path[0] + '/fonts/04B_24__.TTF', 8)
MED_FONT = ImageFont.truetype(TTF, 32)
LRG_FONT = ImageFont.truetype(TTF,  64)
SMALL_FONT_WIDTH, SMALL_FONT_HEIGHT = get_font_dim(SMALL_FONT)
MED_FONT_WIDTH, MED_FONT_HEIGHT = get_font_dim(MED_FONT)

IMG_BG = (255,255,255)

TICK_BG = (0,50,0)
TICK_FG = (255,255,255)
TICK_OUTLINE = TICK_BG
TICK_UNIT_HEIGHT = SMALL_FONT_HEIGHT
RULER_HEIGHT = (3 * TICK_UNIT_HEIGHT) + MED_FONT_HEIGHT

# TILE
BASES_PER_TILE = 32
TILE_SIZE = BASES_PER_TILE * SMALL_FONT_WIDTH

DEFAULT_SEGMENT_WIDTH = BASES_PER_TILE * (2**5)
DEFAULT_IMAGE_WIDTH = TILE_SIZE * (2**5)

## READ PLACEMENT
def sort_reads(reads):
    ''' Sort reads by read start, and determine extremities of contig,
    including read overhangs. '''
    sorted_reads = [(reads[rn].padded_start,
                     reads[rn].padded_end,
                     rn) for rn in reads]
    sorted_reads.sort()
    #    sorted_reads = [r[1] for r in sorted_reads]
    return sorted_reads

def adjust_contig_coordinates(reads):
    ''' adjust \'contig\' start and end to includerelative to read '''
    return (min([reads[rn].padded_start for rn in reads]),
            max([reads[rn].padded_end for rn in reads]))

def stack_reads(reads):
    ''' stack reads ... was tile reads '''
    pass
    
def get_segments(read_start, read_end, segments):
    ''' Get segments that overlap the read '''

    in_segments = []
    for i, segment in enumerate(segments):
        
        segment_start, segment_end, segment_size = segment
        if (segment_start <= read_start <= segment_end) or \
           (segment_start <= read_end <= segment_end) or \
           (segment_start >= read_start and segment_end <= read_end):
            in_segments.append(i)
    return in_segments

def calculate_segments(asm_width):
    ''' Calculates how to segment up the contig. '''
    
    i = 0
    segments = []
    
    if DEFAULT_SEGMENT_WIDTH >= asm_width:
        
        segments.append([0, (asm_width - 1), asm_width])
    else:
        
        while i * DEFAULT_SEGMENT_WIDTH < asm_width:

            segment_start = i * DEFAULT_SEGMENT_WIDTH
            segment_end = segment_start + DEFAULT_SEGMENT_WIDTH
            if segment_end >= asm_width:
                segment_end = asm_width
            segment_size = segment_end - segment_start
            segments.append([segment_start, segment_end, segment_size])
            i += 1

    if segments[-1][2] < BASES_PER_TILE: 
        # Merge small right-most segment and second-to-last segment
        s = segments.pop(-1)
        segments[-1][1] += s[2]
        segments[-1][2] += s[2]
    
    if len(segments) > 1 and (len(segments) % 2) != 0:
        # if odd segment count, merge last and second-to-last
        s = segments.pop(-1)
        segments[-1][1] += s[2]
        segments[-1][2] += s[2]
    
    return segments


def tile_reads(reads, contig):
    ''' Given hash of reads, determines row placement so reads do not overlap.
    Returns array of tuples (read_name, row)'''

    
    INTER_READ_PAD = 5 # space between reads, maybe use to display read name in
                       # future?
    
    current_ends = [0] # stores current end of right-most read in row
    
    read_start, read_end, read_name = reads.pop(0)
    read = contig.reads[read_name]
    current_ends[0] = read_end # + len(read.name)# current end of right-most read in first row
    tiled_reads = [(read, 0)] # stores tiles reads
    while reads: # while we still have reads to process
        read_start, read_end, read_name = reads.pop(0)
        read = contig.reads[read_name]
        is_placed = False
        for (i, row_end) in enumerate(current_ends):
            if read_start > (row_end + INTER_READ_PAD): # if read fits on row
                current_ends[i] = read_end # + len(read.name)# update current end
                tiled_reads.append((read, i)) # append to tiled reads list
                is_placed = True # mark currently processed read as placed
                break
        if not is_placed: # if current read not placed ...
            current_ends.append(read_end) #+ len(read.name)) # add entry to current_ends
            i = len(current_ends) - 1 # how may rows now?
            tiled_reads.append((read, i)) # add read to tiled reads
            
    return tiled_reads, len(current_ends)

def tile_paired_reads(reads, contig):
    ''' Given hash of reads, determines row placement so reads do not overlap.
    Returns array of tuples (read_name, row)'''

    
    INTER_READ_PAD = 5 # space between reads, maybe use to display read name in
                       # future?
                       
    linkers = []
    current_ends = [0] # stores current end of right-most read in row
    read_names = [r[2] for r in reads]
    read_start, read_end, read_name = reads.pop(0)
    read = contig.reads[read_name]
    current_ends[0] = read_end # current end of right-most read in first row
    tiled_reads = [(read, 0)] # stores tiles reads
    if read.is_paired:
        read_names = [r[2] for r in reads]
        paired_read_name = contig.reads[read_name].base_name + '_right'
        if '_right' in read_name:
            paired_read_name = contig.reads[read_name].base_name + '_left'
            
        if contig.reads.get(paired_read_name):
            paired_read = contig.reads[paired_read_name]
            idx = read_names.index(paired_read_name)
            paired_read_start, paired_read_end, paired_read_name = reads.pop(idx)
            current_ends[0] = paired_read_end 
            tiled_reads.append((paired_read, 0))
            linkers.append((read.padded_start, paired_read.padded_end, 0, read.base_name))
    
    while reads: # while we still have reads to process
        read_start, read_end, read_name = reads.pop(0)
        read = contig.reads[read_name]
        is_placed = False
        for (i, row_end) in enumerate(current_ends):
            if read_start > (row_end + INTER_READ_PAD): # if read fits on row
                current_ends[i] = read_end # update current end
                tiled_reads.append((read, i)) # append to tiled reads list
                if read.is_paired:
                    read_names = [r[2] for r in reads]
                    paired_read_name = contig.reads[read_name].base_name + '_right'
                    if '_right' in read_name:
                        paired_read_name = contig.reads[read_name].base_name + '_left'
                    if contig.reads.get(paired_read_name):
                        paired_read = contig.reads[paired_read_name]
                        idx = read_names.index(paired_read_name)
                        paired_read_start, paired_read_end, paired_read_name = reads.pop(idx)
                        current_ends[i] = paired_read_end
                        tiled_reads.append((paired_read, i))
                        linkers.append((read.padded_start, paired_read.padded_end, i, read.base_name))
                
                is_placed = True # mark currently processed read as placed
                break
        if not is_placed: # if current read not placed ...
            current_ends.append(read_end) # add entry to current_ends
            i = len(current_ends) - 1 # how may rows now?
            tiled_reads.append((read, i)) # add read to tiled reads
            if read.is_paired:
                paired_read_name = contig.reads[read_name].base_name + '_right'
                read_names = [r[2] for r in reads]
                if '_right' in read_name:
                    paired_read_name = contig.reads[read_name].base_name + '_left'
                if contig.reads.get(paired_read_name):
                    paired_read = contig.reads[paired_read_name]
                    idx = read_names.index(paired_read_name)
                    paired_read_start, paired_read_end, paired_read_name = reads.pop(idx)
                    current_ends[i] = paired_read_end
                    tiled_reads.append((paired_read, i))
                    linkers.append((read.padded_start, paired_read.padded_end, i, read.base_name))
    return tiled_reads, len(current_ends), linkers


    
def draw_linkers(draw, tiled_reads, contig, y_offset, contig_start, segment_start, segment_size, pair_status):
    ''' draw linkers for paierd end reads '''
    s_x1 = segment_start #* SMALL_FONT_WIDTH
    s_x2 = s_x1 + (segment_size)# * SMALL_FONT_WIDTH)

    for (link_start, link_end, row, template) in tiled_reads:
        bgcolor, fgcolor = assembly.READPAIRSTATUS_TO_RGB[pair_status[template].status]
        ins_size = "%0.1f" % (float(link_end - link_start)/1000)
        y = (row * SMALL_FONT_HEIGHT) + y_offset
        #x1 = ((link_start - contig_start) * SMALL_FONT_WIDTH) - s_x1
        #x2 = ((link_end - contig_start) * SMALL_FONT_WIDTH) - s_x1
        x1 = ((link_start - contig_start)) - s_x1
        x2 = ((link_end - contig_start)) - s_x1
        
        if (x2 > 0) or (x1 < s_x2):
            draw.rectangle([(x1, y), (x2, (y + SMALL_FONT_HEIGHT))], fill = bgcolor)
        for i in range(link_start, link_end):
            if (i % 100) == 0:
                #x = ((i - contig_start) * SMALL_FONT_WIDTH) - s_x1
                x = ((i - contig_start)) - s_x1
                draw.text((x, y), ins_size, fill = fgcolor, font=SMALL_FONT)
            
def draw_linkers_2(draw, tiled_reads, contig, y_offset, contig_start, segment_start, segment_size, pair_status, pxbp):
    ''' draw linkers for paierd end reads '''
    s_x1 = segment_start * pxbp
    s_x2 = s_x1 + (segment_size) * pxbp
    for (link_start, link_end, row, template) in tiled_reads:
        bgcolor, fgcolor = assembly.READPAIRSTATUS_TO_RGB[pair_status[template].status]
        ins_size = "%0.1f" % (float(link_end - link_start)/1000)
        y = (row * SMALL_FONT_HEIGHT) + y_offset
        x1 = ((link_start - contig_start) * pxbp) - s_x1
        x2 = ((link_end - contig_start) * pxbp) - s_x1
        
        if (x2 > 0) or (x1 < s_x2):
            cy = (SMALL_FONT_HEIGHT / 2) + 1
            y1 = y + cy - 1
            y2 = y + cy + 1
            draw.rectangle([(x1, y), (x2, y + SMALL_FONT_HEIGHT - 1)], fill = bgcolor)
            #draw.rectangle([(x1, y), ((x1+4), (y + SMALL_FONT_HEIGHT))], fill = bgcolor)
            #draw.rectangle([((x2-4), y), (x2, (y + SMALL_FONT_HEIGHT))], fill = bgcolor)

            txt = "%s" % (template,)
            tx = x1 + (((x2-x1) - (len(template) * SMALL_FONT_WIDTH)) / 2)
            draw.text((tx, (y+1)), txt, font = SMALL_FONT, fill = (0,0,0))
        #for i in range(link_start, link_end):
        #    if (i % 100) == 0:
        #        x = ((i - contig_start) * pxbp) - s_x1
        #        draw.text((x, y), ins_size, fill = fgcolor, font=SMALL_FONT)



def segment_reads(reads, segments, contig_start):
    ''' Segment reads into groups by chunk of SEGMENT_SIZE.
    Requied to print large images. '''

    segmented_reads = []
    for i in range(len(segments)): segmented_reads.append([])
    
    for read, row in reads:
        read_start = read.padded_start - contig_start
        read_end = read_start + read.num_padded_bases - 1
        in_segments = get_segments(read_start, read_end, segments)
        for i in in_segments: segmented_reads[i].append((read, row))

    return segmented_reads

###########
## REFACTOR CODE BELOW
###########

def draw_consensus(draw, contig, contig_start, segment_size, y_offset):
    ''' Draw consensus sequence '''
    
    # compute coordinates of contig subseq to process bases on segment
    offset = 0
    if contig_start > 0: contig_start -= 1
    if contig_start < 0:
        offset = abs(contig_start) + 1
        contig_start = 0
    contig_end = contig_start + segment_size + 1
    if contig_end > len(contig.consensus): contig_end = len(contig.consensus)

    at_base = len([x for x in contig.consensus[:contig_start] if x != '*'])
    
    # background color for consensus
    draw.rectangle((0, y_offset - 1,segment_size*SMALL_FONT_WIDTH, \
                        y_offset + SMALL_FONT_HEIGHT-1), fill=(255,255,165))
    
    for i in range(contig_start, contig_end):
        n = contig.consensus[i]
        if n != '*': at_base += 1
        x = ((i - contig_start) + offset) * SMALL_FONT_WIDTH
        p = i % 10
        b = at_base % 10
        font_color = (0,0,0)
        if n == '*': font_color = (165,165,165)
        draw.text((x, y_offset), n, font=SMALL_FONT, fill=font_color)


def draw_ruler_2(draw, contig, contig_start, segment_size, y_offset, otype):
    
    # compute coordinates of contig subseq to process based on segment
    offset = 0
    if contig_start > 0: contig_start -= 1
    if contig_start < 0:
        offset = abs(contig_start) + 1
        contig_start = 0
    contig_end = contig_start + segment_size + 1
    if contig_end > len(contig.consensus): contig_end = len(contig.consensus)

    at_base = len([x for x in contig.consensus[:contig_start] if x != '*'])
    
    sfw = SMALL_FONT_WIDTH
    if otype == 'templates':
        sfw = 1
    for i in range(contig_start, contig_end):
        n = contig.consensus[i]
        if n != '*': 
            at_base += 1
            x = ((i - contig_start) + offset) * sfw
            if otype == 'templates':
                if (at_base % 1000 == 0):
                    label = str(at_base)
                    text_width, text_height = draw.textsize(label, font=MED_FONT)
                    draw.text(((x - (text_width / 2)), y_offset), label, font=MED_FONT, fill=(0,0,0))
                    draw.rectangle((x, MED_FONT_HEIGHT,
                                    (x + sfw), RULER_HEIGHT), fill=(0,0,0), outline=IMG_BG)

            if otype == 'templates':
                continue

            if (at_base % 100 == 0):
                label = str(at_base)
                text_width, text_height = draw.textsize(label, font=MED_FONT)
                draw.text(((x - (text_width / 2)), y_offset), label, font=MED_FONT, fill=(0,0,0))
                draw.rectangle((x, MED_FONT_HEIGHT,
                                (x + sfw), RULER_HEIGHT), fill=(0,0,0), outline=IMG_BG)
            elif (at_base % 10) == 0:
                bg = (165,165,165)
                fg = (255,255,255)
                if (at_base % 50) == 0:
                    bg = (0,0,0)
                    
                draw.rectangle(((x - 1), (MED_FONT_HEIGHT + TICK_UNIT_HEIGHT),
                                (x + sfw - 1), RULER_HEIGHT), fill=bg, \
                               outline=bg)
                draw.text((x, (RULER_HEIGHT - (SMALL_FONT_HEIGHT * 2))), str((at_base % 100) / 10), \
                          font=SMALL_FONT, fill=fg)
                draw.text((x, (RULER_HEIGHT - SMALL_FONT_HEIGHT)), "0", font=SMALL_FONT, fill=fg)
            else:
                draw.text((x, (RULER_HEIGHT - SMALL_FONT_HEIGHT)), str(at_base % 10), font=SMALL_FONT, fill=(165,165,165))
    
def draw_ruler_3(draw, contig, contig_start, segment_size, y_offset, otype, pxbp):
    
    # compute coordinates of contig subseq to process based on segment
    offset = 0
    if contig_start > 0: contig_start -= 1
    if contig_start < 0:
        offset = abs(contig_start) + 1
        contig_start = 0
    contig_end = contig_start + segment_size + 1
    if contig_end > len(contig.consensus): contig_end = len(contig.consensus)

    at_base = len([x for x in contig.consensus[:contig_start] if x != '*'])
    
    for i in range(contig_start, contig_end):
        n = contig.consensus[i]
        if n != '*': 
            # at_base += 1
            x = math.ceil(((i - contig_start) + offset) * pxbp)
            if otype == 'templates':
                if (at_base % 1000 == 0):
                    label = str(at_base)
                    text_width, text_height = draw.textsize(label, font=MED_FONT)
                    draw.text(((x - (text_width / 2)), y_offset), label, font=MED_FONT, fill=(0,0,0))
                    draw.rectangle((x, MED_FONT_HEIGHT,
                                    (x + 5), RULER_HEIGHT), fill=(0,0,0), outline=IMG_BG)
                
            at_base += 1

                
def draw_reads(draw, reads, contig_start, segment_start, segment_size, \
                   contig, y_offset, read_status, otype):
    ''' Draw tiled reads '''

    fg = (0,0,0)
    segment_end = segment_start + segment_size - 1
    draw_text = draw.text
    for read, row in reads:
        read_pos = read.padded_start - contig_start
        offset = read_pos - segment_start
        read_start = None
        if offset < 0:
            read_start = abs(offset)
            offset = 0
        else:
            read_start = 0
            
        x1 = offset * SMALL_FONT_WIDTH
        y1 = (row * SMALL_FONT_HEIGHT) + y_offset
        x2 = offset * SMALL_FONT_WIDTH + len(read.sequence[read_start:]) * \
            SMALL_FONT_WIDTH
        y2 = y1  + SMALL_FONT_HEIGHT + 1
        
        try:
            bgcolor, fgcolor = read_status[read.name]
        except KeyError:
            try:
                if read.is_paired:
                    bgcolor, fgcolor = read_status[read.base_name + \
                                                   read.paired_end]
                else:
                    bgcolor, fgcolor = read_status[read.base_name]
            except KeyError:
                bgcolor, fgcolor = assembly.READSTATUS_TO_RGB['Assembled']
                
        if otype == 'readnames':
            read_name_str = ''
            draw.rectangle((x1-1, y1-1, x2, y2-3), fill = bgcolor)
            while (len(read_name_str) + len(read.name) + 5) < read.padded_length():
                read_name_str += (" " * 5) + read.name
            
            read_name_str += (" " * (read.padded_length() - len(read_name_str)))
            draw_text((x1, y1), read_name_str[read_start:], font=SMALL_FONT, fill = fgcolor)
              
        if otype == 'assembly':
            draw.rectangle((x1-1, y1-1, x2, y2-3), fill = (0,0,0))
            draw_text((x1, y1), read.sequence[read_start:], font=SMALL_FONT, fill = (165,165,165))
           
            y = (row * SMALL_FONT_HEIGHT) + y_offset
            cpi = read.padded_start + read_start - 1
            uc = contig.consensus.upper()
            for i, n in enumerate(read.sequence[read_start:]):
                cp = cpi + i
                if cp > 0 and cp < len(contig.consensus):
                    cb = uc[cp]
                    if cb != n:
                        x = (i + offset) * SMALL_FONT_WIDTH
                        draw_text((x, y), n, font=SMALL_FONT, fill = BASE_COLORS[n])

def draw_paired_reads(draw, reads, contig_start, segment_start, segment_size, \
                   contig, y_offset, read_status, pair_status):
    ''' Draw tiled reads '''

    fg = (0,0,0)
    segment_end = segment_start + segment_size - 1
    
    for read, row in reads:
        read_pos = read.padded_start - contig_start
        
        offset = read_pos - segment_start
        read_start = None
        if offset < 0:
            read_start = abs(offset)
            offset = 0
        else:
            read_start = 0
            
        x1 = offset * SMALL_FONT_WIDTH
        y1 = (row * SMALL_FONT_HEIGHT) + y_offset
        x2 = offset * SMALL_FONT_WIDTH + len(read.sequence[read_start:]) * \
            SMALL_FONT_WIDTH
        y2 = y1  + SMALL_FONT_HEIGHT + 1

        bgcolor = (0,0,0)

        try:
            bgcolor = read_status[read.name]
        except KeyError:
            try:
                if read.is_paired:
                    bgcolor = read_status[read.base_name + read.paired_end]
                else:
                    bgcolor = read_status[read.base_name]
            except KeyError:
                bgcolor = (0,0,0)

        bgcolor = (255,255,210)
        if read.is_paired:
            if read.paired_contig == contig.name:
                bgcolor = (0,0,255)
            elif (not read.other_contig):
                bgcolor = (255,0,0)
            else:
                bgcolor = (0,255,0)
        else:
            bgcolor = (255,255,210)
         
        
        fgcolor, bgcolor = assembly.READPAIRSTATUS_TO_RGB[pair_status[read.base_name].status]    
        draw.rectangle((x1-1, y1, x2, y2-3), fill = bgcolor)
        read_name_str = ''
        while (len(read_name_str) + len(read.name) + 5) < read.padded_length():
            read_name_str += (" " * 5) + read.name
        read_name_str += (" " * (read.padded_length() - len(read_name_str)))
        draw.text((x1, y1), read_name_str[read_start:], font=SMALL_FONT, fill = fgcolor)
            
        

        # draw.text((x1, y1), read.sequence[read_start:], font=SMALL_FONT, fill = fgcolor)
        
        # MISMATCHES
        if False:
            for i, n in enumerate(read.sequence[read_start:]):
                x =  (i + offset) * SMALL_FONT_WIDTH
                y = (row * SMALL_FONT_HEIGHT) + y_offset
                cp = read.padded_start + read_start + i - 1
                if cp > 0 and cp < len(contig.consensus):
                    cb = contig.consensus[cp]
                    if cb != n:
                        draw.text((x, y), n, font=SMALL_FONT, fill=(255,0,0))



                    
def get_tick_box_coords(pos, txt, fn, y1):
    
    txt_width, txt_height = ImageDraw.Draw(\
        Image.new("RGB", (1000,1000))).textsize(txt, font=fn)
    print txt_width, txt_height, txt
    y2 = y1 + txt_height
    x1 = pos * SMALL_FONT_WIDTH
    x2 = x1 + txt_width
    return [x1, y1, x2, y2, txt]

def tick_label_coords(contig, contig_start, contig_end):
    
    cfg = ImageConfig()
    pos = 0
    bp = 0
    cfg.rh = 0
    cfg.ticks = {}
    for i, tick_at in enumerate(ticks_at):
        if contig.length() >= tick_at:
            cfg.ticks[tick_at] = {}
            cfg.ticks[tick_at]['font'] = ticks_font[i]
            cfg.ticks[tick_at][''] = ticks_font[i]
            cfg.ticks[tick_at]['ticks'] = []
            cfg.ticks[tick_at]['y'] = cfg.rh
            cfg.rh += get_font_dim(cfg.ticks[tick_at]['font'])[1]
            print "HERE", i, tick_at, cfg.ticks[tick_at]['y'], cfg.rh
            
            
    cfg.by = cfg.rh
    cfg.rh += int(math.log10(contig.length()) + 1) * SMALL_FONT_HEIGHT
    cfg.yc = cfg.rh
    for i in range(contig_start, contig_end):
        
        if i > 0 and i < len(contig.consensus):
            if contig.consensus[i] != '*':
                bp += 1
            for j, tick_at in enumerate(cfg.ticks):
                if (bp % tick_at) == 0:
                    print j
                    label = "%03d Kb" % ((bp / 1000),)
                    tick_box_coords = \
                        get_tick_box_coords(pos, label, \
                                                cfg.ticks[tick_at]['font'],\
                                                cfg.ticks[tick_at]['y'])
                    cfg.ticks[tick_at]['ticks'].append(tick_box_coords)
        pos += 1
        
    return cfg


def draw_ticks(draw, segment, tick_labels):
    segment_start, segment_end, segment_size = segment
    sx1, sx2 = segment_start * SMALL_FONT_WIDTH, segment_end * SMALL_FONT_WIDTH
    for tick_at in tick_labels.ticks:
        tick = tick_labels.ticks[tick_at]
        for (x1, y1, x2, y2, label) in tick['ticks']:
            if (x1 in range(sx1,sx2)) or (x2 in range(sx1, sx2)) or \
                    (sx1 in range(x1, x2) and sx2 in range(x1, x2)):
                print tick_at, sx1, sx2, x1 ,x2, 
                tx1, tx2 = x1 - sx1, x2 - sx1
                
                draw.rectangle((tx1, y1, tx2, y2), fill=TICK_BG, outline=TICK_OUTLINE)
                draw.text((tx1, tick['y']), label, font=tick['font'], \
                              fill=TICK_FG)


def get_shotgun_reads(reads, contig):
    return [r for r in reads]
    # return [r for r in reads if not contig.reads[r[2]].is_paired]

def get_paired_reads(reads, contig):
    pr = []
    for r in reads:
        read = contig.reads[r[2]]
        if read.is_paired:
            paired_read_name = read.base_name + "_right"
            if read.paired_end == '_right':
                paired_read_name = read.base_name + "_left"
            if contig.reads.get(paired_read_name):
                pr.append(r)

    return pr

def draw_segment(img_width, img_height, contig, segment, otype, segmented_shotgun_reads, i, read_status, pair_status, image_name):
    
    segment_start, segment_end, segment_size = segment
    contig_start = contig.start + segment_start
    img = Image.new("RGB", (img_width, img_height), IMG_BG)
    draw = ImageDraw.Draw(img)
    draw_ruler_2(draw, contig, contig_start, segment_size, 0, otype)
        
    if otype != 'templates':
        draw_consensus(draw, contig, contig_start, segment_size, RULER_HEIGHT)
        
    if otype in ('assembly', 'readnames'):
        draw_reads(draw, segmented_shotgun_reads[i], contig.start,
                   segment_start, segment_size, contig,
                   (RULER_HEIGHT + SMALL_FONT_HEIGHT), read_status, otype)
    elif otype == 'templates' and (len(tiled_paired_reads) > 0):
        draw_linkers(draw, linkers, contig,
                     (RULER_HEIGHT + SMALL_FONT_HEIGHT),
                     contig.start,
                     segment_start, segment_size, pair_status)
    img.save(image_name + ".png")



def draw_contig(contig, out_dir, read_status, otype, pair_status, PROCESSES):
    ''' Sort & Tile reads, get image dimensions '''
    
    contig.start, contig.end = adjust_contig_coordinates(contig.reads)
    # stacked_reads = stack_reads(contig.reads)
    sorted_reads = sort_reads(contig.reads)
    asm_width = contig.end - contig.start
    shotgun_reads = get_shotgun_reads(sorted_reads, contig)
    paired_reads = get_paired_reads(sorted_reads, contig)
    # print len(shotgun_reads), len(paired_reads), len(sorted_reads)
    # tiled_reads, asm_height = tile_reads(sorted_reads, contig)
    tiled_shotgun_reads, shotgun_asm_height = tile_reads(sorted_reads, contig)
    paired_asm_height = 0
    tiled_paired_reads, paired_asm_height, linkers = [], 0, None
    if len(paired_reads) > 0:
        tiled_paired_reads, paired_asm_height, linkers = tile_paired_reads(paired_reads, contig)
        
    asm_height = shotgun_asm_height + paired_asm_height
    segments = calculate_segments(asm_width)
    segmented_shotgun_reads = segment_reads(tiled_shotgun_reads, segments, \
                                    contig.start)
    segmented_paired_reads = None
    if len(tiled_paired_reads) > 0:
        segmented_paired_reads = segment_reads(tiled_paired_reads, segments, \
                                               contig.start)
    # tick_labels = tick_label_coords(contig, contig_start, \
    # contig_end)
    pan_width = asm_width * SMALL_FONT_WIDTH
    pan_height = ((asm_height + 1) * SMALL_FONT_HEIGHT) # + tick_labels.rh 
    pan_height += RULER_HEIGHT
    img_height = pan_height
    if otype == 'templates':
        print pxbp
        pan_width = math.ceil(asm_width * pxbp)
        pan_height = ((asm_height + 1) * SMALL_FONT_HEIGHT) # + tick_labels.rh 
        pan_height += RULER_HEIGHT
        img_height = pan_height
    
    images = []
    multiprocessing.freeze_support()
    #pool = multiprocessing.Pool(PROCESSES)
    TASKS = []
    for i, segment in enumerate(segments):
        image_name = out_dir + "%s" % (i)
        images.append(image_name)
        segment_start, segment_end, segment_size = segment
        contig_start = contig.start + segment_start
        img_width = segment_size * SMALL_FONT_WIDTH
        if otype == 'templates':
            img_width = segment_size
        
        #TASKS.append((draw_segment, (img_width, img_height, contig, segment, otype, segmented_shotgun_reads, i, read_status, pair_status, image_name)))
        draw_segment(img_width, img_height, contig, segment, otype, segmented_shotgun_reads, i, read_status, pair_status, image_name)
        # draw_segment(img_width, img_height, contig, segment, otype, segmented_shotgun_reads, i, read_status, pair_status, image_name)
        
        ## img = Image.new("RGB", (img_width, img_height), IMG_BG)
        ## draw = ImageDraw.Draw(img)
        ## draw_ruler_2(draw, contig, contig_start, segment_size, 0, otype)
        
        ## if otype != 'templates':
        ##     draw_consensus(draw, contig, contig_start, segment_size, RULER_HEIGHT)

        ## if otype in ('assembly', 'readnames'):
        ##     draw_reads(draw, segmented_shotgun_reads[i], contig.start, \
        ##                    segment_start, segment_size, contig, \
        ##                    (RULER_HEIGHT + SMALL_FONT_HEIGHT), read_status, otype)
        ## elif otype == 'templates' and (len(tiled_paired_reads) > 0):
        ##     draw_linkers(draw, linkers, contig, \
        ##                  (RULER_HEIGHT + SMALL_FONT_HEIGHT),
        ##                  contig.start, \
        ##                  segment_start, segment_size, pair_status)
        ## img.save(image_name + ".png")
        
    #pool.map_async(calculatestar, TASKS).get(9999999)
    return images, pan_width, pan_height

def draw_contig_2(contig, out_dir, read_status, otype, pair_status, pxbp):
    ''' Sort & Tile reads, get image dimensions '''


    contig.start, contig.end = adjust_contig_coordinates(contig.reads)
    # stacked_reads = stack_reads(contig.reads)
    
    
    sorted_reads = sort_reads(contig.reads)
    asm_width = contig.end - contig.start
    print 'ASM_WIDTH'
    print asm_width
    shotgun_reads = get_shotgun_reads(sorted_reads, contig)
    paired_reads = get_paired_reads(sorted_reads, contig)
    print len(shotgun_reads), len(paired_reads), len(sorted_reads)
    
    
    # tiled_reads, asm_height = tile_reads(sorted_reads, contig)
    tiled_shotgun_reads, shotgun_asm_height = tile_reads(sorted_reads, contig)
    paired_asm_height = 0
    tiled_paired_reads, paired_asm_height, linkers = [], 0, None
    if len(paired_reads) > 0:
        tiled_paired_reads, paired_asm_height, linkers = tile_paired_reads(paired_reads, contig)
        
    asm_height = shotgun_asm_height + paired_asm_height
    
    segments = calculate_segments(asm_width)

    segmented_shotgun_reads = segment_reads(tiled_shotgun_reads, segments, \
                                    contig.start)

    segmented_paired_reads = None
    
    if len(tiled_paired_reads) > 0:
        segmented_paired_reads = segment_reads(tiled_paired_reads, segments, \
                                               contig.start)
    
    # tick_labels = tick_label_coords(contig, contig_start, \
    # contig_end)

    pan_width = asm_width * SMALL_FONT_WIDTH
    pan_height = ((asm_height + 1) * SMALL_FONT_HEIGHT) # + tick_labels.rh 
    pan_height += RULER_HEIGHT
    img_height = pan_height
    
    if otype == 'templates':
        pan_width = math.ceil(asm_width * pxbp)
        pan_height = ((asm_height + 1) * SMALL_FONT_HEIGHT) # + tick_labels.rh 
        pan_height += RULER_HEIGHT
        img_height = pan_height


    images = []

    for i, segment in enumerate(segments):
        
        segment_start, segment_end, segment_size = segment
        contig_start = contig.start + segment_start
        
        img_width = math.ceil(segment_size * pxbp)
        if otype == 'templates':
            img_width = int(math.ceil(segment_size * pxbp))
        print "IMG WIDTH"
        print img_width
        img = Image.new("RGB", (img_width, img_height), IMG_BG)
        
        draw = ImageDraw.Draw(img)
        
        draw_ruler_3(draw, contig, contig_start, segment_size, 0, otype, pxbp)
        
        if otype != 'templates':
            draw_consensus(draw, contig, contig_start, segment_size, RULER_HEIGHT)

        if otype in ('assembly', 'readnames'):
            draw_reads(draw, segmented_shotgun_reads[i], contig.start, \
                           segment_start, segment_size, contig, \
                           (RULER_HEIGHT + SMALL_FONT_HEIGHT), read_status, otype)
        elif otype == 'templates' and (len(tiled_paired_reads) > 0):
            draw_linkers_2(draw, linkers, contig, \
                         (RULER_HEIGHT + SMALL_FONT_HEIGHT),
                         contig.start, \
                         segment_start, segment_size, pair_status, pxbp)

            #draw_paired_reads(draw, segmented_paired_reads[i], contig.start,
            #                  segment_start, segment_size, contig,
            #                  (RULER_HEIGHT + SMALL_FONT_HEIGHT),
            #                  read_status, pair_status)
            
            
        if False: #len(tiled_paired_reads) > 0:
            draw_linkers(draw, linkers, contig, \
                         ((RULER_HEIGHT + SMALL_FONT_HEIGHT) + \
                          (shotgun_asm_height * SMALL_FONT_HEIGHT)),
                         contig.start, \
                         segment_start, segment_size)
        
        
                   
                   
        
        # draw_ruler(draw, contig, contig_start, segment_size, tick_labels)
        # draw_ticks(draw, segment, tick_labels)

        image_name = out_dir + "%s" % (i)
        img.save(image_name + ".png")
        images.append(image_name)
        
    # exit(0)
    return images, pan_width, pan_height



def resize_images(images):

    for image_name in images:
        img = Image.open(image_name + ".png")
        img_width, img_height = img.size
        r_width = int(round(img_width / 2.0))
        r_height = int(round(img_height / 2.0))
        if r_height == 0: r_height = 1
        r_img = img.resize((r_width, r_height), Image.BILINEAR)
        r_img.save(image_name + ".png")

def merge_image_pair(image0, image1):
    img0 = Image.open(image0 + ".png")
    img1 = Image.open(image1 + ".png")
    img0_width, img0_height = img0.size
    img1_width, img1_height = img1.size
    img = Image.new("RGB", ((img0_width + img1_width), \
                            img0_height), IMG_BG)
    img.paste(img0, (0,0))
    img.paste(img1, (img0_width ,0))
    os.remove(image0 + ".png")
    os.remove(image1 + ".png")
    return img

def merge_images(images, out_dir):
    new_images = []
    if len(images) > 1:

       if (len(images)) % 2 == 1:
            img = merge_image_pair(images[-2], images[-1])
            img.save(images[-2] + ".png")
            images.pop(-1)
            
       for i in range(0,len(images),2):
            img = merge_image_pair(images[i], images[i+1])
            image_name =  out_dir + str(i/2)
            img.save(image_name + ".png")
            new_images.append(image_name)
                
    else:
        new_images = images

    return new_images

def calc_levels(dims): return int(math.ceil(math.log(max(dims), 2)))
        
def tile_image(tile_dir, image_name, last_col):
    img = Image.open(image_name + ".png")
    image_width, image_height =  img.size
    num_cols =  int(math.ceil(float(image_width) / TILE_SIZE))
    num_rows = int(math.ceil(float(image_height) / TILE_SIZE))
    row = 0
    while row < num_rows:
        y = row * TILE_SIZE
        y2 = y + TILE_SIZE
        if y2 > image_height: y2 = image_height - 1
        if y == y2: y2 += 1
        col = 0
        while col < num_cols:
            tile_name = "%s_%s" % ((col + last_col), row)
            x = col * TILE_SIZE
            x2 = x + TILE_SIZE
            if x2 > image_width: x2 = image_width - 1
            if x == x2: x2 += 1
            tile_file = open(tile_dir + tile_name + ".png", "wb")
            tile_image = img.crop((x, y, x2, y2))
            tile_image.save(tile_file)
            col += 1
        row += 1

def create_zoomtig(out_dir, num_levels, images, PROCESSES):
    multiprocessing.freeze_support()
    #pool = multiprocessing.Pool(PROCESSES)
    for level in range(num_levels, -1, -1):
        
        
        tile_dir = out_dir + str(level) + "/"
        if not os.path.exists(tile_dir):
            os.mkdir(tile_dir)
        
        last_col = 0
        c = 0
        col_offsets = []
        
        last_col = 0
        for image_name in images:
            img = Image.open(image_name + ".png")
            image_width, image_height =  img.size
            num_cols =  int(math.ceil(float(image_width) / TILE_SIZE))
            col_offsets.append(last_col)
            last_col += num_cols
        TASKS = []
        for idx, image_name in enumerate(images):
            #TASKS.append((tile_image, (tile_dir, image_name, col_offsets[idx])))
            tile_image(tile_dir, image_name, col_offsets[idx])
        #pool.map_async(calculatestar, TASKS).get(9999999)
        #for idx, image_name in enumerate(images):
        #    tile_image(tile_dir, image_name, col_offsets[idx])
        resize_images(images)
        images = merge_images(images, out_dir)
        
    

def mkdir(x):
    if not os.path.exists(x):
        try:
            os.mkdir(x)
        except OSError:
            # print "Cannot create %s" % x
            exit(0)

            
def process_contig(contig, dzi_dir, read_status, otype, pair_status, num_processors):
    ''' Create a zoomtig from a contig '''
    
    out_dir = dzi_dir + "/" + contig.name + "_files/"
    mkdir(out_dir)
    sys.stderr.write("drawing ")
    images, pan_width, pan_height = draw_contig(contig, out_dir, read_status, otype, pair_status, num_processors)
    image_params = {'height': pan_height/float(pan_width),
                    'padded_start': contig.start,
                    'padded_end': contig.end,
                    'ruler_height': (RULER_HEIGHT + SMALL_FONT_HEIGHT)/float(pan_width),
                    'max_depth': max(contig.depth_histogram())
                    }
    
    num_levels = calc_levels((pan_width, pan_height))
    sys.stderr.write("tiling ")
    create_zoomtig(out_dir, num_levels, images, num_processors)
    xml = '<?xml version="1.0" encoding="UTF-8"?>'\
          '<Image TileSize="%s" Overlap="0" Format="png" '\
          'xmlns="http://schemas.microsoft.com/deepzoom/2008">'\
          '<Size Width="%s" Height="%s"/></Image>' % \
          (TILE_SIZE, pan_width, pan_height)
    return xml, image_params
    

def process_contig_2(contig, dzi_dir, read_status, otype, pair_status, pxbp):
    ''' create a zoomtig from a contig (assembly.py).'''
    
    out_dir = dzi_dir + "/" + contig.name + "_files/"
    mkdir(out_dir)
    
    images, pan_width, pan_height = draw_contig_2(contig, out_dir, read_status, otype, pair_status, pxbp)
    image_params = {'height': pan_height/float(pan_width),
                    'padded_start': contig.start,
                    'padded_end': contig.end,
                    'ruler_height': (RULER_HEIGHT + SMALL_FONT_HEIGHT)/float(pan_width),
                    'max_depth': max(contig.depth_histogram())
                    }
    
    num_levels = calc_levels((pan_width, pan_height))
    create_zoomtig(out_dir, num_levels, images)
    xml = '<?xml version="1.0" encoding="UTF-8"?>'\
          '<Image TileSize="%s" Overlap="0" Format="png" '\
          'xmlns="http://schemas.microsoft.com/deepzoom/2008">'\
          '<Size Width="%s" Height="%s"/></Image>' % \
          (TILE_SIZE, pan_width, pan_height)
    return xml, image_params


if __name__ == "__main__":
    pass

