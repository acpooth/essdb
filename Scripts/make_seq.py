#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#------------------------------
# Name:
# Purpose: Functions for the creation of sequences for a networkx tree structure.
#
# @uthor:   acph - dragopoot@gmail.com
#
# Created:
# Copyright:   (c) acph 2014
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
#------------------------------
"""
Functions for the creation of sequences for a networkx tree structure.
The sequences form by backtacing each node from leavs to the root

"""

import networkx as nx

# some version of networkx --- successors and predecessors functions returns iterators
# insted of lists. So, to patch this ... all the aparitions of this functions
# where converted to lists


def seqs_from_leaf(grafod, leaf):
    """
        Regresa las secuencias de pasos passos que van desde un nodo externo u hoja hasta
        el nodo raiz. El grafo debe ser drigido y tener forma de árbol
        Las secuencias estan representadas en forma de lista 
    """
    diam = nx.diameter(grafod.to_undirected()) + 2
    count = 0
    leaf = [leaf]
    seqs = []
    seqs.append(leaf)
    while count < diam:
        for seq in seqs:
            if list(grafod.predecessors(seq[0])) == []:
                continue
            elif len(list(grafod.predecessors(seq[0]))) == 1:
                next = list(grafod.predecessors(seq[0]))
                seq.insert(0, next[0])
            elif len(list(grafod.predecessors(seq[0]))) > 1:
                next = list(grafod.predecessors(seq[0]))
                temp = seq[:]
                for i, n in enumerate(next):
                    if i == 0:
                        seq.insert(0, next[0])
                    else:
                        seqs.append([n] + temp)
        count += 1
    return seqs


def leafs_list(grafod):
    """
        Regresa una lista de los nodos mas externos u hojas.
        Requiere un grafo dirigido.
        Un nodo hoja es aquel que no tienen nodos sucesores.
    """
    leafs = []
    for n in grafod.nodes():
        if list(grafod.successors(n)) == []:
            leafs.append(n)
    return leafs


def all_seq(grafod):
    """
        Regresa una lista con todas las secuencias de pasos desde el nodo de inicio hasta
        cada nodo de termino en un grafo dirigido en forma de árbol.
        Las secuencias estan representadas en forma de lista
    """
    seqs = []
    leafs = leafs_list(grafod)
    for leaf in leafs:
        seqs += seqs_from_leaf(grafod, leaf)
    return seqs


def tfs(grafo, inicio, inicios=[]):
    """
       Tree From Start
       Regresa un árbol (grafo dirigido) a paritir de un nodo de inicio. 
       Toma como argumento opcional una lista de inicios. Si se introduce la lista de inicios,
       la construccion del arbol termina cuando encuentra un nodo de inicio 
    """
    if not grafo.is_directed():
        grafo = grafo.to_directed()
    tree = nx.DiGraph()
    frente = []
    frente.append(inicio)
    while frente != []:
        trans = frente[:]
        frente = []
        for n in trans:
            if len(tree) == 0:
                seqs = [n]
            else:
                seqs = seqs_from_leaf(tree, n)
                seqs = set([s for seq in seqs for s in seq])
            for s in list(grafo.successors(n)):
                if s != n and s not in seqs and s not in inicios:
                    #                if s not in tree: # salida similar a BFS
                    tree.add_edge(n, s)
                    frente.append(s)
    return tree


def seqs_from_starts_bfs(grafo, inicios):
    """
        Regresa una lista con todas las secuencias de pasos desde cada nodo deinicio hasta
        cada nodo de termino en un grafo.
        Usa el algoritmo Breadth First Search (networkx) para construir los arboles a partir de lso cuales se generan
        las secuencias
        Las secuencias estan representadas en forma de lista
    """
    seqs = []
    for i in inicios:
        if grafo.has_node(i):
            tree = nx.traversal.bfs_tree(grafo, i)
            seqs += all_seq(tree)
    return seqs


def seqs_from_starts_tfs(grafo, inicios, ti=True):
    """
        Regresa una lista con todas las secuencias de pasos desde cada nodo deinicio hasta
        cada nodo de termino en un grafo.
        Usa la funcion tfs  (Tree from start, acph) para construir los arboles a partir de los cuales se generan
        las secuencias
        Las secuencias estan representadas en forma de lista
        ti = tomar en cuenta inicios.
    """
    if ti == True:
        lista_inicios = inicios
    else:
        lista_inicios = []
    seqs = []
    for i in inicios:
        if grafo.has_node(i):
            tree = tfs(grafo, i, lista_inicios)
            seqs += all_seq(tree)
    return seqs


if __name__ == '__main__':
    exit()
