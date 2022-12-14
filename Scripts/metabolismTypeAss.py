#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#------------------------------
# Name:     metabolismTypeAss.py
# Purpose:  Updates a ESS database by adding the metabolism type for each sequence.
#
# @uthor:   acph - dragopoot@gmail.com
#
# Created:
# Copyright:   (c) acph 2014
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
#------------------------------

"""
    This script updates the seqs.db adding metabolism type for each sequence

Usage:
     $ python metabolismTypeAss.py [seqs.db]
"""

import sqlite3 as s3

#sqlite3 utilities
def updateDB(db, table, change, where=''):
    """
        Updates a sqlite database, no commit. You need to commit after al your updates be done.
    """
    db.execute("""UPDATE %s SET %s %s""" % (table, change, where))
    return db

def add_column(db, table, columndef):
    """
        Takes a sqlite3 db object and adds a defined column to the selected table.
        columndef = sqlite3 column definition (name datatype)
    """
    db.execute("""ALTER TABLE %s ADD COLUMN %s """ % (table, columndef))
    return db

#END sqlite3 utilities    


def updateMetabolismType(db, keggmaps):
    """
        Update seqs.db adding Metabolism type to the metabolism field using the keggmaps dictionary that. keggmaps dictionary keys are the metabolic mapas id; and the value of each item is a list wich first value is the metabolic map name and the second  value is the metabolic category
    """
    for k, v in keggmaps.iteritems():
        updateDB(db, 'seqs', 'metabolism="%s"' % v[1], 'WHERE mapid LIKE "___%s"' % k)
        updateDB(db, 'map_stat', 'metabolism="%s"' % v[1], 'WHERE mapid LIKE "___%s"' % k)
    db.commit()
    return db # not necesary return


keggmaps = {'00720': ['Carbon fixation pathways in prokaryotes', 'Energy Metabolism'],
            '00410': ['beta-Alanine metabolism', 'Metabolism of Other Amino Acids'],
            '00680': ['Methane metabolism', 'Energy Metabolism'],
            '00920': ['Sulfur metabolism', 'Energy Metabolism'],
            '00603': ['Glycosphingolipid biosynthesis - globo series', 'Glycan Biosynthesis and Metabolism'],
            '00601': ['Glycosphingolipid biosynthesis - lacto and neolacto series', 'Glycan Biosynthesis and Metabolism'],
            '00600': ['Sphingolipid metabolism', 'Lipid Metabolism'],
            '00604': ['Glycosphingolipid biosynthesis - ganglio series', 'Glycan Biosynthesis and Metabolism'],
            '00550': ['Peptidoglycan biosynthesis', 'Glycan Biosynthesis and Metabolism'],
            '00311': ['Penicillin and cephalosporin biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '00310': ['Lysine degradation', 'Amino Acid Metabolism'],
            '00312': ['beta-Lactam resistance', 'Biosynthesis of Other Secondary Metabolites'],
            '00010': ['Glycolysis / Gluconeogenesis', 'Carbohydrate Metabolism'],
            '00250': ['Alanine, aspartate and glutamate metabolism', 'Amino Acid Metabolism'],
            '00253': ['Tetracycline biosynthesis', 'Metabolism of Terpenoids and Polyketides'],
            '00540': ['Lipopolysaccharide biosynthesis', 'Glycan Biosynthesis and Metabolism'],
            '00910': ['Nitrogen metabolism', 'Energy Metabolism'],
            '00140': ['Steroid hormone biosynthesis', 'Lipid Metabolism'],
            '00750': ['Vitamin B6 metabolism', 'Metabolism of Cofactors and Vitamins'],
            '00300': ['Lysine biosynthesis', 'Amino Acid Metabolism'],
            '00020': ['Citrate cycle (TCA cycle)', 'Carbohydrate Metabolism'],
            '00380': ['Tryptophan metabolism', 'Amino Acid Metabolism'],
            '00642': ['Ethylbenzene degradation', 'Biodegradation and Metabolism'],
            '00904': ['Diterpenoid biosynthesis', 'Metabolism of Terpenoids and Polyketides'],
            '00905': ['Brassinosteroid biosynthesis', 'Metabolism of Terpenoids and Polyketides'],
            '00906': ['Carotenoid biosynthesis', 'Metabolism of Terpenoids and Polyketides'],
            '00900': ['Terpenoid backbone biosynthesis', 'Metabolism of Terpenoids and Polyketides'],
            '00901': ['Indole alkaloid biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '00902': ['Monoterpenoid biosynthesis', 'Metabolism of Terpenoids and Polyketides'],
            '00903': ['Limonene and pinene degradation', 'Metabolism of Terpenoids and Polyketides'],
            '00908': ['Zeatin biosynthesis', 'Metabolism of Terpenoids and Polyketides'],
            '00909': ['Sesquiterpenoid biosynthesis', 'Metabolism of Terpenoids and Polyketides'],
            '00860': ['Porphyrin and chlorophyll metabolism', 'Metabolism of Cofactors and Vitamins'],
            '00531': ['Glycosaminoglycan degradation', 'Glycan Biosynthesis and Metabolism'],
            '00533': ['Glycosaminoglycan biosynthesis - keratan sulfate', 'Glycan Biosynthesis and Metabolism'],
            '00532': ['Glycosaminoglycan biosynthesis - chondroitin sulfate', 'Glycan Biosynthesis and Metabolism'],
            '00534': ['Glycosaminoglycan biosynthesis - heparan sulfate', 'Glycan Biosynthesis and Metabolism'],
            '00430': ['Taurine and hypotaurine metabolism', 'Metabolism of Other Amino Acids'],
            '00740': ['Riboflavin metabolism', 'Metabolism of Cofactors and Vitamins'],
            '00660': ['C5-Branched dibasic acid metabolism', 'Carbohydrate Metabolism'],
            '00331': ['Clavulanic acid biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '00330': ['Arginine and proline metabolism', 'Amino Acid Metabolism'],
            '00980': ['Metabolism of xenobiotics by cytochrome P450', 'Biodegradation and Metabolism'],
            '00981': ['Insect hormone biosynthesis', 'Metabolism of Terpenoids and Polyketides'],
            '00982': ['Drug metabolism - cytochrome P450', 'Biodegradation and Metabolism'],
            '00983': ['Drug metabolism - other enzymes', 'Biodegradation and Metabolism'],
            '00130': ['Ubiquinone and other terpenoid-quinone biosynthesis', 'Metabolism of Cofactors and Vitamins'],
            '00232': ['Caffeine metabolism', 'Biosynthesis of Other Secondary Metabolites'],
            '00230': ['Purine metabolism', 'Nucleotide Metabolism'],
            '00231': ['Puromycin biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '00030': ['Pentose phosphate pathway', 'Carbohydrate Metabolism'],
            '00626': ['Naphthalene degradation', 'Biodegradation and Metabolism'],
            '00522': ['Biosynthesis of 12-, 14- and 16-membered macrolides', 'Metabolism of Terpenoids and Polyketides'],
            '00523': ['Polyketide sugar unit biosynthesis', 'Metabolism of Terpenoids and Polyketides'],
            '00520': ['Amino sugar and nucleotide sugar metabolism', 'Carbohydrate Metabolism'],
            '00521': ['Streptomycin biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '00524': ['Butirosin and neomycin biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '00281': ['Geraniol degradation', 'Metabolism of Terpenoids and Polyketides'],
            '00280': ['Valine, leucine and isoleucine degradation', 'Amino Acid Metabolism'],
            '00770': ['Pantothenate and CoA biosynthesis', 'Metabolism of Cofactors and Vitamins'],
            '00670': ['One carbon pool by folate', 'Metabolism of Cofactors and Vitamins'],
            '00440': ['Phosphonate and phosphinate metabolism', 'Metabolism of Other Amino Acids'],
            '00120': ['Primary bile acid biosynthesis', 'Lipid Metabolism'],
            '00121': ['Secondary bile acid biosynthesis', 'Lipid Metabolism'],
            '00040': ['Pentose and glucuronate interconversions', 'Carbohydrate Metabolism'],
            '00960': ['Tropane, piperidine and pyridine alkaloid biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '00966': ['Glucosinolate biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '00965': ['Betalain biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '00351': ['DDT degradation', 'Biodegradation and Metabolism'],
            '00350': ['Tyrosine metabolism', 'Amino Acid Metabolism'],
            '00643': ['Styrene degradation', 'Biodegradation and Metabolism'],
            '00195': ['Photosynthesis', 'Energy Metabolism'],
            '00640': ['Propanoate metabolism', 'Carbohydrate Metabolism'],
            '00290': ['Valine, leucine and isoleucine biosynthesis', 'Amino Acid Metabolism'],
            '00514': ['Other types of O-glycan biosynthesis', 'Glycan Biosynthesis and Metabolism'],
            '00513': ['Various types of N-glycan biosynthesis', 'Glycan Biosynthesis and Metabolism'],
            '00512': ['Mucin type O-Glycan biosynthesis', 'Glycan Biosynthesis and Metabolism'],
            '00511': ['Other glycan degradation', 'Glycan Biosynthesis and Metabolism'],
            '00510': ['N-Glycan biosynthesis', 'Glycan Biosynthesis and Metabolism'],
            '00592': ['alpha-Linolenic acid metabolism', 'Lipid Metabolism'],
            '00591': ['Linoleic acid metabolism', 'Lipid Metabolism'],
            '00590': ['Arachidonic acid metabolism', 'Lipid Metabolism'],
            '00190': ['Oxidative phosphorylation', 'Energy Metabolism'],
            '00760': ['Nicotinate and nicotinamide metabolism', 'Metabolism of Cofactors and Vitamins'],
            '00196': ['Photosynthesis - antenna proteins', 'Energy Metabolism'],
            '00450': ['Selenocompound metabolism', 'Metabolism of Other Amino Acids'],
            '00052': ['Galactose metabolism', 'Carbohydrate Metabolism'],
            '00053': ['Ascorbate and aldarate metabolism', 'Carbohydrate Metabolism'],
            '00051': ['Fructose and mannose metabolism', 'Carbohydrate Metabolism'],
            '00791': ['Atrazine degradation', 'Biodegradation and Metabolism'],
            '00790': ['Folate biosynthesis', 'Metabolism of Cofactors and Vitamins'],
            '00950': ['Isoquinoline alkaloid biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '00650': ['Butanoate metabolism', 'Carbohydrate Metabolism'],
            '00340': ['Histidine metabolism', 'Amino Acid Metabolism'],
            '01052': ['Type I polyketide structures', 'Metabolism of Terpenoids and Polyketides'],
            '01055': ['Biosynthesis of vancomycin group antibiotics', 'Metabolism of Terpenoids and Polyketides'],
            '01054': ['Nonribosomal peptide structures', 'Metabolism of Terpenoids and Polyketides'],
            '01057': ['Biosynthesis of type II polyketide products', 'Metabolism of Terpenoids and Polyketides'],
            '01056': ['Biosynthesis of type II polyketide backbone', 'Metabolism of Terpenoids and Polyketides'],
            '00500': ['Starch and sucrose metabolism', 'Carbohydrate Metabolism'],
            '01058': ['Acridone alkaloid biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '00062': ['Fatty acid elongation', 'Lipid Metabolism'],
            '00061': ['Fatty acid biosynthesis', 'Lipid Metabolism'],
            '00710': ['Carbon fixation in photosynthetic organisms', 'Energy Metabolism'],
            '00460': ['Cyanoamino acid metabolism', 'Metabolism of Other Amino Acids'],
            '00260': ['Glycine, serine and threonine metabolism', 'Amino Acid Metabolism'],
            '00100': ['Steroid biosynthesis', 'Lipid Metabolism'],
            '00780': ['Biotin metabolism', 'Metabolism of Cofactors and Vitamins'],
            '00785': ['Lipoic acid metabolism', 'Metabolism of Cofactors and Vitamins'],
            '00940': ['Phenylpropanoid biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '00941': ['Flavonoid biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '00942': ['Anthocyanin biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '00943': ['Isoflavonoid biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '00944': ['Flavone and flavonol biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '00945': ['Stilbenoid, diarylheptanoid and gingerol biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '01040': ['Biosynthesis of unsaturated fatty acids', 'Lipid Metabolism'],
            '00471': ['D-Glutamine and D-glutamate metabolism', 'Metabolism of Other Amino Acids'],
            '00472': ['D-Arginine and D-ornithine metabolism', 'Metabolism of Other Amino Acids'],
            '00473': ['D-Alanine metabolism', 'Metabolism of Other Amino Acids'],
            '00625': ['Chloroalkane and chloroalkene degradation', 'Biodegradation and Metabolism'],
            '00071': ['Fatty acid metabolism', 'Lipid Metabolism'],
            '00072': ['Synthesis and degradation of ketone bodies', 'Lipid Metabolism'],
            '00073': ['Cutin, suberine and wax biosynthesis', 'Lipid Metabolism'],
            '00621': ['Dioxin degradation', 'Biodegradation and Metabolism'],
            '00620': ['Pyruvate metabolism', 'Carbohydrate Metabolism'],
            '00623': ['Toluene degradation', 'Biodegradation and Metabolism'],
            '00622': ['Xylene degradation', 'Biodegradation and Metabolism'],
            '00270': ['Cysteine and methionine metabolism', 'Amino Acid Metabolism'],
            '01053': ['Biosynthesis of siderophore group nonribosomal peptides', 'Metabolism of Terpenoids and Polyketides'],
            '00830': ['Retinol metabolism', 'Metabolism of Cofactors and Vitamins'],
            '00624': ['Polycyclic aromatic hydrocarbon degradation', 'Biodegradation and Metabolism'],
            '00627': ['Aminobenzoate degradation', 'Biodegradation and Metabolism'],
            '00401': ['Novobiocin biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '00400': ['Phenylalanine, tyrosine and tryptophan biosynthesis', 'Amino Acid Metabolism'],
            '00402': ['Benzoxazinoid biosynthesis', 'Biosynthesis of Other Secondary Metabolites'],
            '00730': ['Thiamine metabolism', 'Metabolism of Cofactors and Vitamins'],
            '00360': ['Phenylalanine metabolism', 'Amino Acid Metabolism'],
            '00361': ['Chlorocyclohexane and chlorobenzene degradation', 'Biodegradation and Metabolism'],
            '00362': ['Benzoate degradation', 'Biodegradation and Metabolism'],
            '00363': ['Bisphenol degradation', 'Biodegradation and Metabolism'],
            '00364': ['Fluorobenzoate degradation', 'Biodegradation and Metabolism'],
            '00930': ['Caprolactam degradation', 'Biodegradation and Metabolism'],
            '01051': ['Biosynthesis of ansamycins', 'Metabolism of Terpenoids and Polyketides'],
            '00633': ['Nitrotoluene degradation', 'Biodegradation and Metabolism'],
            '00630': ['Glyoxylate and dicarboxylate metabolism', 'Carbohydrate Metabolism'],
            '00564': ['Glycerophospholipid metabolism', 'Lipid Metabolism'],
            '00565': ['Ether lipid metabolism', 'Lipid Metabolism'],
            '00562': ['Inositol phosphate metabolism', 'Carbohydrate Metabolism'],
            '00563': ['Glycosylphosphatidylinositol(GPI)-anchor biosynthesis', 'Glycan Biosynthesis and Metabolism'],
            '00561': ['Glycerolipid metabolism', 'Lipid Metabolism'],
            '00240': ['Pyrimidine metabolism', 'Nucleotide Metabolism'],
            '00480': ['Glutathione metabolism', 'Metabolism of Other Amino Acids'],
            '00970' :  ['Aminoacyl-tRNA biosynthesis', 'Genetic Information Processing; Translation']}



if __name__ == '__main__':
    """
    $ python metabolismTypeAss.py [seqs.db]
    """
    from sys import argv

    print('[INFO] Opening database')
    db = s3.connect(argv[1])
    print('[INFO] Creating Metabolism column')
    add_column(db, "map_stat", "metabolism TEXT")
    print('[INFO] Filling metabolism column')
    updateMetabolismType(db, keggmaps)
    db.close()
    print('[END] Process complete :P')
