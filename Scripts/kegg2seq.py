#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#------------------------------
# Name:     kegg2seq.py
# Purpose:   Creates a database of enzymatic step sequences (ESS).
#
# @uthor:   acph - dragopoot@gmail.com
#
# Created:
# Copyright:   (c) acph 2014
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
#------------------------------

"""
Creates a database of Enzymatic Step Sequences (ESS).

Usage:

   $  python kegg2seq.py


"""
import argparse
from parse_kgml import *
from make_seq import *
from traducir import *
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
import numpy as np
import sqlite3 as s3

from xml.parsers.expat import ExpatError


def argumentParser():
    parser = argparse.ArgumentParser(description='Script for creation of Enzymatic Step Sequence database',
                                     epilog="Please follow the instructions in README.md")
    parser.add_argument('-g', '--graphtype', default='nodirected', help="""Type of graph to use
    for the contruction of ESS. 'directed' or 'nodirected' [nodirected]. Using 'directed'
    graph takes into account the directionality of the reactions for the construction of BFS trees
    and affect the total number of ESS generated.""")
    return parser.parse_args()


def map_stats(react_dict, pathway, grafog, inicios, gen2ec):
    """
        Creates a dictionary whith some stats form a given map.
        Needs to introduce a parser dict, a parsed xml tree, a graph of genes
        a start nodes list and a dictionary gen > EC numbers.
        All of this can be generated whith modules here.

        >>>>>>>>>
        SUGGESTED MODIFICATION : input the kgml file and generate all the
        resources whitin the function.
        Else, This function may remaind same if i create
        a Map Class which can have as methos all the functions in these
        modules
        <<<<<<<<<
    """
    title = pathway.getAttribute('title')
    name = pathway.getAttribute('name').split(':')[1]
    reacciones = len(react_dict)
    reversibles = 0
    genes = []
    genes_sin_ec = []
    EC3 = []
    EC4 = []
    map_links = []
    for k, v in react_dict.iteritems():
        if v['type'] == 'reversible':
            reversibles += 1
        genes += v['gene'].split()
        if v.has_key('maplink'):
            map_links += v['maplink']
    for g in genes:
        if gen2ec.has_key(g):
            E = gen2ec[g]
            EC4 += E
            for e in E:
                e = e.split('.')
                EC3.append('.'.join(E[:3]))
        else:
            genes_sin_ec.append(g)
    resultado = {}
    resultado['Title'] = title
    resultado['Name'] = name
    resultado['Reactions'] = reacciones
    resultado['Reversibles'] = reversibles
    resultado['Genes'] = len(set(genes))
    resultado['EC3'] = len(set(EC3))
    resultado['EC4'] = len(set(EC4))
    resultado['Sin_EC'] = (len(set(genes_sin_ec)), genes_sin_ec)
    resultado['Nodos'] = len(grafog)
    resultado['Inicios'] = len(inicios)
    resultado['Links_mapas'] = len(set(map_links))
    resultado['Links'] = ' '.join(set(map_links))
    return resultado


def save_hist(values_list, map_name, path):
    if values_list == []:
        return 0
    maxi = max(values_list)
    plt.hist(values_list, bins=maxi, range=(
        1, maxi + 1), align='left', label=map_name)
    plt.legend()
    plt.xlabel('Seq lenght')
    plt.ylabel('Count')
    plt.xticks(np.arange(maxi + 1))
    plt.savefig('%s.png' % (path + map_name))
    plt.close()


def createDB(name):
    """
        This function creates a database of gene and EC numbers sequences.
        Table seqs has the following structure:

        id = unique identification number for each genes chain (Primary Key)
                gen_seq = gene sequence. Some positions in the sequence contain 2 or more genes (holoenzymes, multiple protein subunits), gene positions in the sequence are separated by ">"
                ec3 = EC numbers sequence. 3 levels of clasification. Each position separated by ":"
                ec4 = EC numbers sequence. 4 levels of clasification. Each position separated by ":"
                arch = architecture sequence --pendient--
                len = sequence lenght
                map = KEGG metabolic map from which came the gene sequence
                mapid = KEGG id of the metabolic map (Foreign Key refered to Field mapid in sp_maps table)
                metabolism = KEGG metabolic category --Done, another script fills this field--
                sp = species from witch came the sequence. KEGG 3 letter code

       Table map_stat has the folowing structure:
       mapid = KEGG map unique identificator (primary key)
       name = map name
       sp = specie
       reactions = number of reactions
       reversibles = number of reversible reactions
       genes= number of genes
       ec3s = number of ECs 3 levels
       ec4s = number of ECs 4 levels
       no_EC = number of proteins whith no EC assigned
       nodes = number of nodes of the map graph
       starts = number ofstarts nodes from which the sequences were constructed
       maplinks = links to other maps
    """
    db = s3.connect(name)
    db.execute("""CREATE TABLE seqs
(id INTEGER PRIMARY KEY,
gen_seq TEXT,
ec3 TEXT,
ec4 TEXT,
arch TEXT,
len INTEGER,
map TEXT,
mapid TEXT,
metabolism TEXT,
sp TEXT,
FOREIGN KEY (mapid) REFERENCES map_stat (mapid))""")
    db.execute("""PRAGMA foreign_keys = TRUE""")
    db.execute("""CREATE TABLE map_stat
(mapid TEXT PRIMARY KEY,
name TEXT,
sp TEXT,
reactions INTEGER,
reversibles INTEGER,
genes INTEGER,
ec3s INTEGER,
ec4s INTEGER,
no_EC INTEGER,
nodes INTEGER,
starts INTEGER,
maplinks TEXT,
nseqs INTEGER) """)

    return db


if __name__ == '__main__':
    """
    This script nedds to be in kegg2seq/Scripts
    """
    import os
    args = argumentParser()
    # OLD
    # sps = [d for d in os.listdir('../Enzymes/') if '.' not in d]
    # sps = [d for d in sps if '_' not in d]
    # END OLD
    sps = [d for d in os.listdir('../Enzymes/')]
    sps = [d.split('.')[0] for d in sps]
    # removing sps for 'issues'
    # sps.remove('aaa')
    # end removing
    sps.sort()
    # if 'bsu' in sps: # estr if y su contenido se pueden comentar
    #     sps.remove('bsu') # para pasar a bsu al final de la lista
    #     sps.append('bsu')
    currentpath = os.path.realpath('.')
    if not os.path.exists('../Db/'):
        os.mkdir('../Db/')
    # if not os.path.exists('../Output/'):
    #     os.mkdir('../Output/')
    # if not os.path.exists('../Stats/'):
    #     os.mkdir('../Stats/')
    # textdb = open('../Db/seqs.txt', 'w')
    db = createDB('../Db/seqs.db')
    os.chdir('../')
    seqid = 1
    all_len = []
    # count for intermedate commiting
    count = 0
    for sp in sps:
        os.chdir('Maps/%s/' % sp)
        # To create text database
        # os.mkdir('/home/acph/warp/Mis_doc/essdb/Output/%s' % sp)
        # os.mkdir('/home/acph/warp/Mis_doc/essdb/Output/%s/Seqs' % sp)
        # os.mkdir('/home/acph/warp/Mis_doc/essdb/Output/%s/Hist' % sp)

        #######################################
        # WARNING: This line was set manually #
        #######################################
        # g2e = gen2ec('/home/acph/Desktop/essdb/Enzymes/%s.list' % sp )
        g2e = gen2ec('../../Enzymes/%s.list' % sp)
        sp_len = []
        limit = "{}01000".format(sp)
        maps = [m for m in os.listdir('./') if '.xml' in m and m < limit]
        # remove this line
        # maps = [m for m in os.listdir('./') if m < limit]
        maps.sort()

        for m in maps:
            print 'Processing %s map %s' % (sp, m)
            # map_file = open('/home/acph/warp/Mis_doc/essdb/Output/%s/Seqs/%s.dat' % (sp, m[:-4]), 'w')
            try:
                d, d_ = kgml2dict(m)
            except ExpatError:
                with open('kgml2dict.log', 'a') as k2dlog:
                    k2dlog.write(
                        '[WARN] %s map %s ommited with ExpatError\n' % (sp, m))
                    continue
            if args.graphtype == 'directed':
                g = nx.DiGraph(lista_adyacencia_genes(d))
            elif args.graphtype == 'nodirected':
                g = nx.Graph(lista_adyacencia_genes(d))
            else:
                print "Wrong option for graphtype, please remove ../DB folder and try again."
                print "Exiting!"
                exit()

            i = iniciosMG_gene(d, g)
            if i == []:
                i = iniciosMG_gene(d, g, x=5, y=4)
                if i == []:
                    i = iniciosMG_gene(d, g, x=6, y=5)
                    if i == []:
                        i = iniciosMG_gene(d, g, x=7, y=6)
                        if i == []:
                            i = iniciosMG_gene(d, g, x=8, y=7)
            stat = map_stats(d, d_, g, i, g2e)
            seqs = seqs_from_starts_bfs(g, i)
            lenghts = []
            for se in seqs:
                lenghts.append(len(se))
            sp_len += lenghts
            # save_hist(lenghts, m[:-4], '/home/acph/warp/Mis_doc/essdb/Output/%s/Hist/' % sp)
            if lenghts == []:
                Max = 0
                Min = 0
            else:
                Max = max(lenghts)
                Min = min(lenghts)
#             map_file.write(
#             """# Title:\t%s
# # Map:\t%s
# # Reactions:\t %d
# # Reversibles:\t%d
# # Genes:\t%d
# # ECs-3:\t%d
# # ECs-4:\t%d
# # No EC:\t%d
# # Nodes:\t%d
# # Starts:\t%d\t%s
# # MapLinks:\t%d\t%s
# # Sequences :\t%d
# # Max lenght:\t%d
# # Min lenght:\t%d


# """ % (stat['Title'], stat['Name'], stat['Reactions'], stat['Reversibles'], stat['Genes'],stat['EC3'],stat['EC4'], stat['Sin_EC'][0], stat['Nodos'], len(i), ' '.join(i), stat['Links_mapas'], stat['Links'], len(seqs), Max, Min ))
            mid = 1
            db.execute("""INSERT INTO map_stat VALUES ('%s', '%s', '%s', %d, %d, %d, %d, %d, %d, %d, %d, '%s', %d) """ % (
                stat['Name'], stat['Title'], sp, stat['Reactions'],
                stat['Reversibles'], stat['Genes'], stat['EC3'],
                stat['EC4'], stat['Sin_EC'][0], stat['Nodos'], len(i),
                stat['Links_mapas'], len(seqs)))
            print "   ... {} sequences".format(len(seqs))
            for s in seqs:
                gen = '>'.join(s)
                e3 = traducir_seq3(s, g2e)
                e4 = traducir_seq4(s, g2e)
                # map_file.write("%d\t%s\t%s\t%s\n" % (mid, gen, e3, e4))
                # textdb.write("%d\t%s\t%s\t%s\t%d\t%s\t%s\t%s\n" % (seqid, gen, e3, e4, len(s), stat['Title'], stat['Name'], sp  ))
                db.execute(
                    """INSERT INTO seqs (gen_seq, ec3, ec4, len, map, mapid, sp)
VALUES ("%s", "%s", "%s", %i, "%s", "%s", "%s")""" % (gen, e3, e4, len(s), stat['Title'], stat['Name'], sp))

                mid += 1
                seqid += 1
            # map_file.close()
            count += 1
            if count % 2000 == 0:
                print 'Intermedial commit!'
                db.commit()

        all_len += sp_len
        ###########################################################
        # WARNING: This lines also must be specified manually.... #
        # Anoying past me!!!                                      #
        ###########################################################
        # saving organism histogram
        # save_hist(sp_len, sp, '/home/acph/warp/Mis_doc/essdb/Output/%s/' % sp)
        # save_hist(sp_len, sp, '/home/acph/Desktop/essdb/Output/')

        # Return to root
        # os.chdir('/home/acph/Desktop/essdb/')
        os.chdir('../../')

    print 'Commiting'
    db.commit()
    print 'Creating index'
    db.execute("CREATE INDEX id_index ON seqs (id)")
    db.close()
    # textdb.close()
    # os.chdir(currentpath)
    # save_hist(all_len, 'Seqs_DB', '../Stats/')
    print 'Ending!!!! '
