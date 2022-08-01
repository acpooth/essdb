# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 19:52:50 2011

@author: acph

Code for parse kgml and create graphs
some graphs like eco00190 (oxidative phosphorylation) can't  be proccesed using this code.
"""

from xml.dom import minidom
import networkx as nx

# ------------------ parse kgml ----------------------


def kgml2dict(kgml):
    """
        Genera un diccionario a partir de un archivo kgml. El diccionario tiene como llaves el nombre de la reacción y como valores
n        varios atributos de la reaccion, como genes, reactivos y productos
    """
    # leer xml con minidom
    handler = open(kgml, 'r')
    lista = []
    for line in handler.readlines():
        lista.append(line.strip())
    handler.close()
    string = ''
    for l in lista:
        if string == '':
            string += l
        elif l[0] == '<':
            if string[-1] == '>':
                string += l
            else:
                string += ' '
                string += l
        else:
            string += ' '
            string += l
    parsed = minidom.parseString(string)
    path = parsed.childNodes[2]
    # dictionary de reacciones
    react_dict = {}
    for reaction in path.getElementsByTagName('reaction'):
        # if not react_dict.has_key(reaction.getAttribute('name')):
        if reaction.getAttribute('name') not in react_dict:
            react_dict[reaction.getAttribute('name')] = {}
            entry = react_dict[reaction.getAttribute('name')]
            entry['substrate'] = [comp.getAttribute(
                'name') for comp in reaction.getElementsByTagName('substrate')]
            entry['product'] = [comp.getAttribute(
                'name') for comp in reaction.getElementsByTagName('product')]
            entry['type'] = reaction.getAttribute('type')
            entry['id'] = reaction.getAttribute('id')
    # diccionario de traduccion id - rn
    id_rn = {}
    # for k, v in react_dict.items():
    for k, v in react_dict.iteritems():
        id_rn[v['id']] = k
    # crea un diccionario de links a otros mapas y agrega los nombres de los genes a cada reaccion
    map_dict = {}
    for entry in path.getElementsByTagName('entry'):
        # seleciona unicamente los links a mapas de la especie
        if entry.getAttribute('type') == 'map' and 'path:map' not in entry.getAttribute('name'):
            map_dict[entry.getAttribute('id')] = entry.getAttribute('name')
        if entry.getAttribute('type') == 'gene':
            if entry.hasAttribute('reaction'):
                # If Entry not  a reaction ( not in reaction nodes)
                # then skip
                if entry.getAttribute('reaction') not in react_dict:
                    with open('k2s.log', 'a') as outf:
                        line = "kgml2dict---------------\n" +\
                            "Problem with reaction-entry" +\
                            " relation:\n{} \n".format(kgml)
                        outf.write(line)
                    continue
                react_dict[entry.getAttribute(
                    'reaction')]['gene'] = entry.getAttribute('name')
    # agrega los links a otros mapas en la reaccion correspondiente en react_dict
    for link in path.getElementsByTagName('relation'):
        if link.getAttribute('type') == 'maplink':
            if link.getAttribute('entry1') in id_rn:
                if link.getAttribute('entry2') in map_dict:
                    idrn = link.getAttribute('entry1')
                    if "maplink" in react_dict[id_rn[idrn]]:
                        react_dict[id_rn[idrn]]["maplink"].append(
                            map_dict[link.getAttribute('entry2')])
                    else:
                        react_dict[id_rn[idrn]]["maplink"] = [
                            map_dict[link.getAttribute('entry2')]]
            if link.getAttribute('entry2') in id_rn:
                if link.getAttribute('entry1') in map_dict:
                    idrn = link.getAttribute('entry2')
                    if "maplink" in react_dict[id_rn[idrn]]:
                        react_dict[id_rn[idrn]]["maplink"].append(
                            map_dict[link.getAttribute('entry1')])
                    else:
                        react_dict[id_rn[idrn]]["maplink"] = [
                            map_dict[link.getAttribute('entry1')]]
    return react_dict, path


# ------------------ graph ----------------------
def lista_adyacencia(react_dict):
    """
        Crea una lista de adyacencia de reacciones a partir del diccionario creado en kgml2dict()
    """
    lista = []
    for r1k, r1v in react_dict.iteritems():
        if r1v['type'] == "reversible":
            products = r1v['product'] + r1v['substrate']
        else:
            products = r1v['product']
        for p in products:
            for r2k, r2v in react_dict.iteritems():
                if r2v['type'] == 'reversible':
                    substrates = r2v['substrate'] + r2v['product']
                else:
                    substrates = r2v['substrate']
                for s in substrates:
                    if p == s:
                        lista.append((r1k, r2k))
    return lista

# --- genes


def lista_adyacencia_genes(react_dict):
    """
        Crea una lista de adyacencia de genes a partir del diccionario creado en kgml2dict()
    """
    lista = []
    for r1k, r1v in react_dict.iteritems():
        if r1v['type'] == "reversible":
            products = r1v['product'] + r1v['substrate']
        else:
            products = r1v['product']
        for p in products:
            for r2k, r2v in react_dict.iteritems():
                if r2v['type'] == 'reversible':
                    substrates = r2v['substrate'] + r2v['product']
                else:
                    substrates = r2v['substrate']
                for s in substrates:
                    if p == s:
                        lista.append((r1v['gene'], r2v['gene']))
    return lista


def grafoRX(reac_dict):
    g = nx.Graph(lista_adyacencia(reac_dict))
    return g


# ------------------ start nodes list ----------------------
#    substrates = set([s for reac in reac_dict.itervalues() for s in reac['substrate']])

def inicios(reac_dict):
    """
    Crea una lista de reacciones de inicio a partir de un diccionario de reacciones creado por kgml2dict().
    Toma en cuenta la reaversibilidad de las reacciones.
    Inicio = reacción sin producto en el mapa
    """
#    products = set([p for reac in reac_dict.itervalues() for p in reac['product']])
    products = []
    for reac in reac_dict.itervalues():
        for p in reac['product']:
            products.append(p)
        if reac['type'] == 'reversible':
            for s in reac['substrate']:
                products.append(s)
    products = set(products)
    inicios = []
    for rek, rev in reac_dict.iteritems():
        if rev['type'] == 'reversible':
            substrates = rev['substrate'] + rev['product']
        else:
            substrates = rev['substrate']
        for s in substrates:
            if s not in products:
                inicios.append(rek)
    return list(set(inicios))


def iniciosM(reac_dict):
    """
    Crea una lista de reacciones de inicio a partir de un diccionario de reacciones creado por kgml2dict().
    Toma en cuenta la reaversibilidad de las reacciones y las reacciones que tienen vinculos a otros mapas
    Inicio = reacción sin producto en el mapa o reaccion vinculada a otro mapa
    """
    products = []
    for reac in reac_dict.itervalues():
        for p in reac['product']:
            products.append(p)
        if reac['type'] == 'reversible':
            for s in reac['substrate']:
                products.append(s)
    products = set(products)
#    products = set([p for reac in reac_dict.itervalues() for p in reac['product']])
    inicios = []
    for rek, rev in reac_dict.iteritems():
        if rev['type'] == 'reversible':
            substrates = rev['substrate'] + rev['product']
        else:
            substrates = rev['substrate']
        for s in substrates:
            if s not in products:
                inicios.append(rek)
        if 'maplink' in rev:
            # if rev.has_key('maplink'):
            inicios.append(rek)
    return list(set(inicios))


def iniciosMG(reac_dict, grafo, x=4, y=3):
    """
    Crea una lista de reacciones de inicio a partir de un diccionario de reacciones creado por kgml2dict().
    Toma en cuenta la reaversibilidad de las reacciones y las reacciones que tienen vinculos a otros mapas
    Inicio = reacción sin producto en el mapa o reaccion vinculada con otro mapa con una conectividad menor a x para reacciones reversibles y menor a y para reacciones no reversibles
    Default: x = 4, y = 3

    """
    products = []
    for reac in reac_dict.itervalues():
        for p in reac['product']:
            products.append(p)
        if reac['type'] == 'reversible':
            for s in reac['substrate']:
                products.append(s)
    products = set(products)
#    products = set([p for reac in reac_dict.itervalues() for p in reac['product']])
    inicios = []
    for rek, rev in reac_dict.iteritems():
        if rev['type'] == 'reversible':
            substrates = rev['substrate'] + rev['product']
        else:
            substrates = rev['substrate']
        for s in substrates:
            if s not in products:
                inicios.append(rek)
        if 'maplink' in rev:
            # if rev.has_key('maplink'):
            if rev['type'] == 'reversible':
                if len(list(grafo.neighbors(rek))) < x:
                    inicios.append(rek)
            else:
                if len(list(grafo.neighbors(rek))) < y:
                    inicios.append(rek)
    return list(set(inicios))


# ---genes

def inicios_genes(reac_dict):
    """
    Crea una lista de reacciones de inicio a partir de un diccionario de reacciones creado por kgml2dict().
    Toma en cuenta la reaversibilidad de las reacciones.
    Inicio = reacción sin producto en el mapa
    """
   # products = set([p for reac in reac_dict.itervalues() for p in reac['product']])
    products = []
    for reac in reac_dict.itervalues():
        for p in reac['product']:
            products.append(p)
        if reac['type'] == 'reversible':
            for s in reac['substrate']:
                products.append(s)
    products = set(products)
    inicios = []
    for rek, rev in reac_dict.iteritems():
        if rev['type'] == 'reversible':
            substrates = rev['substrate'] + rev['product']
        else:
            substrates = rev['substrate']
        for s in substrates:
            if s not in products:
                inicios.append(rev['gene'])
    return list(set(inicios))


def iniciosM_genes(reac_dict):
    """
    Crea una lista de reacciones de inicio a partir de un diccionario de reacciones creado por kgml2dict().
    Toma en cuenta la reaversibilidad de las reacciones y las reacciones que tienen vinculos a otros mapas
    Inicio = reacción sin producto en el mapa o reaccion vinculada a otro mapa
    """
    products = []
    for reac in reac_dict.itervalues():
        for p in reac['product']:
            products.append(p)
        if reac['type'] == 'reversible':
            for s in reac['substrate']:
                products.append(s)
    products = set(products)
#    products = set([p for reac in reac_dict.itervalues() for p in reac['product']])
    inicios = []
    for rek, rev in reac_dict.iteritems():
        if rev['type'] == 'reversible':
            substrates = rev['substrate'] + rev['product']
        else:
            substrates = rev['substrate']
        for s in substrates:
            if s not in products:
                inicios.append(rev['gene'])
        if 'maplink' in rev:
            # if rev.has_key('maplink'):
            inicios.append(rev['gene'])
    return list(set(inicios))


def iniciosMG_gene(reac_dict, grafo, x=4, y=3):
    """
    Crea una lista de reacciones de inicio a partir de un diccionario de reacciones creado por kgml2dict().
    Toma en cuenta la reaversibilidad de las reacciones y las reacciones que tienen vinculos a otros mapas. Requiere un grafo de genes
    Inicio = reacción sin producto en el mapa o reaccion vinculada con otro mapa con una conectividad menor a x para reacciones reversibles y menor a y para reacciones no reversibles
    Default: x = 4, y = 3

    """
    products = []
    for reac in reac_dict.itervalues():
        for p in reac['product']:
            products.append(p)
        if reac['type'] == 'reversible':
            for s in reac['substrate']:
                products.append(s)
    products = set(products)
#    products = set([p for reac in reac_dict.itervalues() for p in reac['product']])
    inicios = []
    for rek, rev in reac_dict.iteritems():
        if rev['type'] == 'reversible':
            substrates = rev['substrate'] + rev['product']
        else:
            substrates = rev['substrate']
        for s in substrates:
            if s not in products:
                inicios.append(rev['gene'])
        if 'maplink' in rev:
            # if rev.has_key('maplink'):
            if rev['type'] == 'reversible':
                if grafo.has_node(rev['gene']):
                    if len(list(grafo.neighbors(rev['gene']))) < x:
                        inicios.append(rev['gene'])
            else:
                if grafo.has_node(rev['gene']):
                    if len(list(grafo.neighbors(rev['gene']))) < y:
                        inicios.append(rev['gene'])
    return list(set(inicios))


# utiles -------------------------------


"""
for n in g: # quita los nodos que estan conectados solo con si mismos
    if len(g.neighbors(n)) == 1 and  g.neighbors(n)[0] == n:
        g.remove_node(n)
"""

"""
glu, glu_ = kgml2dict('eco00010.xml')
aa , aa_ = kgml2dict('eco00250.xml')
ac , ac_ = kgml2dict('eco00071.xml')
k, k_ = kgml2dict('eco00020.xml')
p, p_ = kgml2dict('eco00230.xml')
py, py_ = kgml2dict('eco00240.xml')

gglu = nx.Graph(lista_adyacencia(glu))
gaa = nx.Graph(lista_adyacencia(aa))
gac = nx.Graph(lista_adyacencia(ac))
gk = nx.Graph(lista_adyacencia(k))
gp = nx.Graph(lista_adyacencia(p))
gpy = nx.Graph(lista_adyacencia(py))

gglug = nx.Graph(lista_adyacencia_genes(glu))
gaag = nx.Graph(lista_adyacencia_genes(aa))
gacg = nx.Graph(lista_adyacencia_genes(ac))
gkg = nx.Graph(lista_adyacencia_genes(k))
gpg = nx.Graph(lista_adyacencia_genes(p))
gpyg = nx.Graph(lista_adyacencia_genes(py))

iglu = iniciosMG(glu, gglu)
iaa = iniciosMG(aa, gaa)
iac = iniciosMG(ac, gac)
ik = iniciosMG(k, gk)
ip = iniciosMG(p, gp)
ipy = iniciosMG(py, gpy)

iglug = iniciosMG_gene(glu, gglug)
iaag = iniciosMG_gene(aa, gaag)
iacg = iniciosMG_gene(ac, gacg)
ikg = iniciosMG_gene(k, gkg)
ipg = iniciosMG_gene(p, gpg)
ipyg = iniciosMG_gene(py, gpyg)

"""

if __name__ == '__main__':
    sys.exit(main())
