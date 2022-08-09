# -*- coding: utf-8 -*-
#
#------------------------------
# Name: traducir.py
# Purpose: Translate a sequence of genes to a sequences of EC numbers (ESS) using
# the first 3 or 4 levels of EC classification.
#
# @uthor:   acph - dragopoot@gmail.com
#
# Created:
# Copyright:   (c) acph 2014
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
#------------------------------

"""
Functions for the translation of a sequence of genes to a sequences of EC numbers (ESS)
using the first 3 or 4 levels of EC classification.
"""

def traducir_seq4(seq , dic) :
    """Translate a sequence of genes to a sequences of EC numbers (4 levels). 
    Manages multiple genes per position in seq.

    seq = string of genes separated by '>'
    dict = translation dictionary of genes to EC

    Return.
    A string of EC numbers separated by colon ':' """
    if type(seq) == type(' '):
        seq = seq.split('>')
    traduccion = []
    for gen in seq :
        if dic.has_key(gen):
            traduccion.append(dic[gen][0])
        elif len(gen.split()) > 1:
            genes = gen.split()
            genes = [g for g in genes if dic.has_key(g)]
            ecs = [ec for gene in genes for ec in dic[gene]]
            if ecs == []:
                e = '9.9.9.9'
            else:
                e = ecs[0]
                for ec in ecs:
                    if ecs.count(ec) == len(genes):
                        e = ec
            traduccion.append(e)
        else:
            traduccion.append('9.9.9.9')
    trad = ':'.join(traduccion)
#    trad += '\n'
    return trad

def traducir_seq3(seq, dic):
    """Translate a sequence of genes to a sequences of EC numbers (3 levels). 
    Manages multiple genes per position in seq.

    seq = string of genes separated by '>'
    dict = translation dictionary of genes to EC

    Return.
    A string of EC numbers separated by colon ':' """
    if type(seq) == type(' '):
        seq = seq.split('>')
    traduccion = []
    for gen in seq :
        if dic.has_key(gen):
            EC = dic[gen][0]
            EC = EC.split('.')
            traduccion.append('.'.join(EC[:3]))
        elif len(gen.split()) > 1:
            genes = gen.split()
            genes = [g for g in genes if dic.has_key(g)]
            ecs = [ec for gene in genes for ec in dic[gene]]
            if ecs == []:
                e = '9.9.9.9'
            else:
                e = ecs[0]
                for ec in ecs:
                    if ecs.count(ec) == len(genes):
                        e = ec
            EC = e.split('.')
            traduccion.append('.'.join(EC[:3]))
        else:
            traduccion.append('9.9.9')
            
    trad = ':'.join(traduccion)
#    trad += '\n'
    return trad

def gen2ec(archivo): #archivo == sp_enzime.list
    """Creates a dictinary for the translation of genes to
    EC number using the kegg sp.list files

    archivo = filename , str"""
    dic = {}
    infile = open(archivo, 'r')
    line = infile.readline()
    while line != '' :
        item = line.split()
        gen = item[0]
        ec = item[1]
        if not dic.has_key(gen):
            dic[gen] = [ec[3:]]
        else:
            dic[gen].append(ec[3:])
        line = infile.readline()
    infile.close()
    return dic

def traducir(archivo, salida3, salida4, dic):
    """Translate the gene sequences in file 'archivo' and stores in EC
    sequences of 3 levels (salida3) and 4 leveles (salida4)  """
    infile = open(archivo, 'r')
    out3 = open(salida3, 'w')
    out4 = open(salida4, 'w')
    linea = infile.readline()
    out3.write(linea)
    out4.write(linea)
    while linea != '':
        if linea[0].islower():
            out3.write(traducir_seq3(linea,dic))
            out4.write(traducir_seq4(linea,dic))
        linea = infile.readline()
    infile.close()
    out3.close()
    out4.close()
        
if __name__ == '__main__':
    pass
    # dic = gen2ec('../Enzimas/eco/eco_enzyme.list')
    # traducir('./eco00030_caminos.txt', '30_3EC.txt', '30_4EC.txt', dic)
