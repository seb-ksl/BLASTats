#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re
import urllib.request
from Bio.Blast import NCBIWWW, NCBIXML


def fetch_organism(string):
    """
    Returns a list containing organisms names contained in BLAST match title.
    Names are filtered to ignore names < 2 words (like "Bacillus cereus"), names containing "group" and names
    repeated several times in a single BLAST match title.
    """

    pattern = re.compile('\[(.*?)\]')
    matchlist = []
    for match in pattern.findall(string):
        if match not in matchlist and "group" not in match and len(match.split()) > 2:
            matchlist.append(match)

    # organism = matchlist[0]
    # for match in matchlist:
    #     if len(match.split()) > 2:
    #         if match.split()[2] != "group":
    #             organism = match
    #             break

    return matchlist


def fetch_genomes_quantity():
    """
    Returns a dictionary containing the number of sequenced genomes available at NCBI for Bc, Bt, Ba and Bs.
    """
    try:

        genomes_quantity = {}
        for organism in ['bacillus+cereus', 'bacillus+thuringiensis', 'bacillus+anthracis', 'bacillus+subtilis']:
            ans_bc = urllib.request.urlopen("http://www.ncbi.nlm.nih.gov/genome/?term={}".format(organism))
            page_bc = ans_bc.read()
            # Because ans.read() returns a byte string while regex findall takes unicode strings
            page_bc = page_bc.decode()

            # What is retrieved from NCBI's page here are 3 occurences of numbers between brackets being respectively
            # the number of genomes, scaffolds/contigs and SRA/traces
            pattern = re.compile('\[(\d+)\]')
            grep = pattern.findall(page_bc)
            # Convert strings list to ints list
            for i in range(len(grep)):
                grep[i] = int(grep[i])

            if organism == 'bacillus+cereus':
                genomes_quantity['bc'] = sum(grep)
            elif organism == 'bacillus+thuringiensis':
                genomes_quantity['bt'] = sum(grep)
            elif organism == 'bacillus+anthracis':
                genomes_quantity['ba'] = sum(grep)
            elif organism == 'bacillus+subtilis':
                genomes_quantity['bs'] = sum(grep)

        return genomes_quantity

    except urllib.request.URLError:
        print("Error:\nCould not reach NCBI's website.\nPlease check your internet connection and retry.")
        sys.exit(0)


def blast():
    """
    BLASTs protein sequence at NCBI and writes the result in a "blast_results.xml" file in the execution directory.
    """
    seq = sys.argv[sys.argv.index('--www') + 1]
    query = NCBIWWW.qblast('blastp', 'nr', seq, hitlist_size=500)
    of_ = open('blast_results.xml', 'w')
    of_.write(query.read())
    of_.close()
    query.close()


def main():
    verbose = False
    if '-v' in sys.argv:
        verbose = True

    if '--www' in sys.argv:
        if verbose:
            print("BLASTing sequence at NCBI...")
        blast()

    # Lists containing hits for each member
    total, bacilli, bc, bt, ba, bs, sp = [], [], [], [], [], [], []
    # Dictionary containing number of sequenced genomes for each member
    if verbose:
        print("Retrieving genomes quantities at NCBI...")
    genomes = fetch_genomes_quantity()
    identity_threshold = .75

    if verbose:
        print("Parsing BLAST results...")
    results_handle = open('blast_results.xml')
    results = NCBIXML.read(results_handle)
    for result in results.alignments:
        for hsp in result.hsps:
            if hsp.identities / results.query_letters > identity_threshold:
                organism_list = fetch_organism(result.title)
                for organism in organism_list:
                    if organism not in total:
                        total.append(organism)
                        if "Bacillus" in organism:
                            bacilli.append(organism)
                            if "cereus" in organism:
                                bc.append(organism)
                            elif "thuringiensis" in organism:
                                bt.append(organism)
                            elif "anthracis" in organism:
                                ba.append(organism)
                            elif "subtilis" in organism:
                                bs.append(organism)
                            else:
                                sp.append(organism)
                break

    print('==============================')
    print("Total: {}".format(len(total)))
    print("|-Bacilli: {}".format(len(bacilli)))
    print("|--- cereus: {0} ({1}%)".format(len(bc), int(len(bc) / genomes['bc'] * 100)))
    print("|--- thuringiensis: {0} ({1}%)".format(len(bt), int(len(bt) / genomes['bt'] * 100)))
    print("|--- anthracis: {0} ({1}%)".format(len(ba), int(len(ba) / genomes['ba'] * 100)))
    print("|--- subtilis: {0} ({1}%)".format(len(bs), int(len(bs) / genomes['bs'] * 100)))
    print("|--- Other Bacilli: {}".format(len(sp)))
    print("|-Other non-Bacilli: {}".format(len(total)-len(bacilli)))
    print('==============================')

    if '-l' in sys.argv:
        bc.sort()
        bt.sort()
        ba.sort()
        bs.sort()
        sp.sort()
        print("\nBacillus cereus:\n-----------------------")
        for bac in bc:
            print(bac)
        print("\nBacillus thuringiensis:\n-----------------------")
        for bac in bt:
            print(bac)
        print("\nBacillus anthracis:\n-----------------------")
        for bac in ba:
            print(bac)
        print("\nBacillus subtilis:\n-----------------------")
        for bac in bs:
            print(bac)
        print("\nOther Bacilli:\n-----------------------")
        for bac in sp:
            print(bac)


if __name__ == '__main__':
    main()
