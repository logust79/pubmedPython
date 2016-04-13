#!/usr/bin/env python
'''
Scrape the table from RetNet disease, and store the data in JSON
'''
import re
import urllib2
import json
from optparse import OptionParser
from bs4 import BeautifulSoup

###############################
# parse options
###############################
usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-u", "--url",
                  dest="url", default='https://sph.uth.edu/Retnet/disease.htm',
                  help="Retnet url: [default: %default]")
parser.add_option("-o", "--output",
                  dest="output", default='../retnet.json',
                  help="output file [default: %default]")

(options, args) = parser.parse_args()
##############################
##############################

# get the url
page = urllib2.urlopen(options.url).read()
soup = BeautifulSoup(page, "html.parser")

# parse
retnet = {}
for tr in soup.find_all('tr')[2:]:
    tds = tr.find_all('td')
    
    # skip all irrelevant rows
    if len(tds) < 3 or not tds[0].a.b:
        continue
    
    # tds[0] encodes genes, and omim ids/links
    # tds[2] encodes diseases
    gene_name = tds[0].a.b.string

    # some genes don't have omim ids
    a_ = tds[0].find_all('a')
    omim = []
    if len(a_) > 1:
        # has omim id
        for a in a_[1:]:
            omim.append(a.string)

    # disease description
    disease = tds[2].get_text().replace('[Gene]','').strip()
    # get mode, u:unknown, x:X-linked, d:dominant, r:recessive, m:mitochondrion, dr: dominant or recessive
    mode = 'u'
    if re.search('recessive', disease):
        if re.search('dominant', disease):
            mode = 'dr'
        else:
            mode = 'r'
    elif re.search('dominant', disease):
        mode = 'd'
    else:
        # check location
        loc = tds[1].string
        if loc:
            if loc[0] == 'X':
                mode = 'x'
            elif loc[:4] == 'mito':
                mode = 'm'
    # write to retnet
    if gene_name in retnet:
        raise 'repetitive gene_name: %s?' % gene_name
    retnet[gene_name] = {'mode': mode, 'omim': omim, 'disease': disease}

# write to output
outf = open(options.output, 'w')
print len(retnet)
json.dump(retnet, outf)
