#!/usr/bin/python/ 
import re 
import argparse
import cairo
import numpy as np

#classes: Sequence, Motif, Artist, Parser

parser = argparse.ArgumentParser()

parser.add_argument("-f", help="fasta filename")
parser.add_argument("-m", help="motif filename")

args = parser.parse_args()

IUPAC = {"u": "T", "n":"ACTG", "y":"CT", "r": "AG","b":"CTG", "d":"AGT", "h":"ACT", "v":"ACG", "w":"AT", "s":"CG", "m":"AC", "k":"GT"}
palette = {1:[90,214,255],2:[253,130,138],3:[190,255,105],4:[182,137,255],5:[92,255,171]}
linecol = [0,0,0]
rectcol = [0,0,0]
motiflocs = {}

def convertrgb(rgb_list):
            '''this function converts rgb 0-255 to rgb 0-1 values'''
            rgb = np.array(rgb_list)/255
            return rgb

def locatemotifs(header, sequence, dimotifs):
    loclist = []
    for motif in dimotifs:
        restring = dimotifs[motif].upper()
        locs = [loc.start() for loc in re.finditer(restring, sequence.upper())]
        if locs != []:
            loclist.append((motif,locs))
    motiflocs[header] = loclist
    return motiflocs

class Artist:
    def __init__(self, di):
        '''
        This class represents an artist where they will draw a line of length "total"
        and a rectangle of length "exon"
        di[header] = (seqlen, exonloc, motiflocs)
        '''
        ## Data ##
        self.name = args.f.split('.')[0]        
        self.di = di
        self.start = 0
    ## Methods ##
    def test(self):   
        print(self.di)

    def draw(self):
        #setting up total dimensions info
        dim_info = []
        for hdr in self.di:
            dim_info.append(self.di[hdr][0])
        WIDTH = max(dim_info) + 200
        HEIGHT = len(self.di) * 80

        #establish canvas
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
        ctx = cairo.Context(surface)

        #line
        ctx.set_line_width(1)

        # Set a background color
        ctx.save()
        ctx.set_source_rgb(1,1,1)
        ctx.paint()
        ctx.restore()

        motifcolors = []
        for i, motif in enumerate(di_motifs):
            motifcolors.append((motif, convertrgb(palette[i+1])))
            

        for i, hdr in enumerate(self.di):
            i = i+1
            moveto = i * 70
            ctx.set_source_rgba(linecol[0], linecol[1], linecol[2], 1)
            ctx.move_to(100,moveto)
            ctx.line_to(100+self.di[hdr][0],moveto)
            ctx.stroke()

            #rectangle - exon
            ctx.rectangle(100+self.di[hdr][1][0], moveto - 10, self.di[hdr][1][1], 20)  # Rectangle(x0, y0, x1, y1)
            ctx.set_source_rgb(rectcol[0], rectcol[1], rectcol[2])
            ctx.fill()

            #rectangles - motifs
            for i, motiftup in enumerate(self.di[hdr][2]):
                motifstarts = motiftup[1]
                motifname = motiftup[0]
                motiflen = len(motiftup[0])

                for tup in motifcolors:
                    if tup[0] == motifname:
                        for start in motifstarts:
                            ctx.rectangle(100+start, moveto - 10, motiflen, 20)  # Rectangle(x0, y0, x1, y1)
                            ctx.set_source_rgba(tup[1][0], tup[1][1], tup[1][2], .70)
                            ctx.fill()        
        #legend
        for i, name in enumerate(motifcolors):
            ctx.set_source_rgb(name[1][0], name[1][1], name[1][2])
            ctx.move_to(10 + (i * 50),10)
            ctx.show_text(name[0].upper())
        for i,hdr in enumerate(self.di):
            ctx.set_source_rgba(linecol[0], linecol[1], linecol[2], 1)
            ctx.move_to(10, (50+ i * 70))
            ctx.show_text(hdr)

        surface.write_to_png(f"{self.name}.png")  # Output to PNG

class Sequence:
    def __init__(self, seqstring, header):
        '''This class represents a sequence where a string gets passed and has
        functions that provide its introns and exons attributes by
        parsing through the string'''
        ## Data ##
        self.string = seqstring
        self.header = header
        self._introns = []
        self._exon = []
        self.strlen = len(self.string)

    ## Methods ##
    def return_hdr(self):
        return self.header

    def parse_seq(self):
        intron = ''
        exon = ''
        for nt in self.string:
            #first intron capture
            if nt.islower() == True and exon == '':
                intron = intron + nt
            #first nt of exon, add first intron to attribute self._introns
            elif nt.isupper() == True and exon == '':
                exon = exon + nt
                self._introns.append(intron)
                intron = ''
            #add rest of exon
            elif nt.isupper() == True:
                exon = exon + nt
            #intron after exon, add exon to attribute self._exon
            else:
                self._exon = exon
                intron = intron + nt
        self._introns.append(intron)

        exonpos = (len(self._introns[0]), len(self._exon))
        return (self.strlen, exonpos)

class Motif:
    def __init__(self, motifstring):
        '''This class represents a sequence where a string gets passed and has
        functions that provide its introns and exons attributes by
        parsing through the string'''
        ## Data ##
        self.string = motifstring.lower()
        self.restring = ''
        self.motiflist = []
        
    ## Methods ##
    def motif_restring(self):
        for nt in self.string:
            if nt in IUPAC:
                self.restring = self.restring + f"[{IUPAC[nt]}]"
            else:
                self.restring = self.restring + nt
        return self.restring

        
class Parser:
  def __init__(self, mfile, ffile):
    '''This class represents the parsing of the fasta file with sequences and motif
    file with a list of motifs and can output a collection of png maps of each
    of the sequences'''

    ## Data ##
    self.ffile = ffile
    self.mfile = mfile

  ## Methods ##
  def parse_mfile(self):
    with open(self.mfile, 'r') as mfh:
        mlist = []
        allmotifs = {}
        for line in mfh:
            line = line.strip()
            mlist.append(line)
        for motif in mlist:
            allmotifs[motif] = Motif(motif).motif_restring()
    return allmotifs
    
  def parse_ffile(self):
    with open(self.ffile, 'r') as ffh:
        flist = []
        hlist = []
        seq = ''
        header = ''

        for line in ffh:
            if ">" in line:
                line = line.strip()
                header = line
                hlist.append(header)
                flist.append(seq)
                seq = ''
                continue
            line = line.strip()
            seq = seq + line
        flist.append(seq)
    hdr_seq_di = dict(zip(hlist, flist[1:]))
    return hdr_seq_di

#Parsing motif file and fasta file 
parsing = Parser(args.m, args.f)

#Initializing Motif Objects with Regex string of Motif as value and Motifs as key
di_motifs = parsing.parse_mfile()

#Initializing Sequence Objects and Populating a Dictionary of Sequence Info
di_hdr_seq = parsing.parse_ffile()
di_seq_info = {}
for hdr in di_hdr_seq:
    di_seq_info[hdr] = Sequence(di_hdr_seq[hdr], hdr).parse_seq()

#Locating Motif Start Sites in Sequence
for hdr in di_hdr_seq:
    locatemotifs(hdr, di_hdr_seq[hdr], di_motifs)

#Feeding Appropriate Sequence Info to Artist Class for Figure Making
artist_di = {}
for i, hdr in enumerate(motiflocs):
    artist_di[hdr] = (di_seq_info[hdr][0], di_seq_info[hdr][1], motiflocs[hdr])
Artist(artist_di).draw()



