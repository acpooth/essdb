# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 17:03:45 2011

@author: acph

Traduce de seq de genes a secuencias de EC numbers de 3 y 4 niveles
NO QUITA LOS redundantes aun
"""

def traducir_seq4(seq , dic) :
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
    dic = gen2ec('../Enzimas/eco/eco_enzyme.list')
    traducir('./eco00030_caminos.txt', '30_3EC.txt', '30_4EC.txt', dic)
