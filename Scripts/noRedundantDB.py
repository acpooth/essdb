#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#------------------------------
# Name:     noRedundantDB.py
# Purpose:   Creates a table of non redundat ESS in a ESS database.
#
# @uthor:   acph - dragopoot@gmail.com
#
# Created:
# Copyright:   (c) acph 2014
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
#------------------------------

"""
This script creates a new table in the ESS vanilla database
named nrseqs that contains the non redundant EC3 ESS (nrESS).
This table contains the nrid of each nrESS and the EC3 ESS.
Additionally, creates a column nrid in the table seqs. This
column is vinculated with the nrid column in nrseqs.

NOTE: the main difference with the legacy script is that
this script creates a more simple database. allowing the
correspondance of the important data from each nrESS
to be obtained in seqs table and eliminating the redundant
data created in legacy.

IMPORTANT. To create nr table, its necesary to previously update
the database with the metabolism type assigned to each sequece,
using metabolismTypeAss.py script

Usage 
   $ python noRedundantDB.py seqs.db
"""

import sqlite3 as s3

# sqlite utilities


def distinct_elements(db, column, table, where=''):
    """Uses a SQLite DB objet and creates a tuple of unique elements of
    field 'column' of a 'table'.
    Optional parameter where = sqlite where statement for refine the qwery"""
    query = db.execute("""SELECT DISTINCT %s FROM %s %s""" %
                       (column, table, where))
    unique = [element[0] for element in query.fetchall()]
    unique.sort()
    return tuple(unique)

# sqlite utilities - END


def create_nr_table(db, table_name='nrseqs', ec='ec3'):
    """Creates a table nrseqs of non redundant sequences in seqs.db of ec3 or ec4.

    Arguments:
    - `db`:
    - `table_name`:
    - `ec`: [ec3]|ec4
    """
    ######################################################################
    # NOTE                                                               #
    # Maybe ommit the len value, this value can easily be computed after #
    ######################################################################
    db.execute("""CREATE TABLE %s
               (nrid INTEGER PRIMARY KEY,
               %s TEXT UNIQUE,
               len INTEGER)""" % (table_name, ec))
    db.execute("""CREATE INDEX nrid_index ON %s (nrid)""" % table_name)
    db.commit()
    return db


def create_nrid_column(db):
    """Creates a nrid column in seqs table. The column contains the
    corresponding nrid from nrseqs

    Arguments:
    - `db`: ESS database sqlite3 database object
    """
    columns_seqs = db.execute('PRAGMA table_info(seqs)')
    columns_seqs = columns_seqs.fetchall()
    columns_seqs = [i[1] for i in columns_seqs]
    if not 'nrid' in columns_seqs:
        db.execute("""ALTER TABLE seqs ADD COLUMN nrid INTEGER""")
    return db


def fill_nr(db, table_name='nrseqs', ec='ec3'):
    """
        This function fills the no redundant table.
    """
    nr = distinct_elements(db, ec, 'seqs')
    lnr = len(nr)  # nr number of sequences
    count = 1

    #######################################################################
    # NOTES                                                               #
    #                                                                     #
    # 1. Maybe create a new table that contains only the correspondence   #
    #    between id and nrid                                              #
    # 2. The seqs id is indexed, so maybe, it is easier (and faster) to   #
    #    assign the correspondence of nrid in ess tableby indexing with   #
    #    te id instead? but requiere some kind of for loop to assign with #
    #    each id , mmmm                                                   #
    # 3. the other may may be to create a index for ec3 in seqs...        #
    #    sounds good. However, how many extra memmory this will need      #
    #######################################################################

    for seq in nr:
        print "Prosessing sequence {0}/{1}...".format(count, lnr)
        lenght = len(seq.split(':'))
        try:
            db.execute("INSERT INTO {} ({}, len) VALUES ('{}', {})".format(
                table_name, ec, seq, lenght))
        except:
            print table_name, ec, seq, lenght
            exit()
        db.execute("UPDATE seqs SET nrid={} WHERE ec3='{}'".format(count, seq))
        count += 1
    db.commit()
    return db  # this return is not necesary


if __name__ == '__main__':
    """$ noRedundatnDB.py [seqs.db]
"""
    from sys import argv

    db = s3.connect(argv[1])
    create_nr_table(db, 'nrseqs')
    create_nrid_column(db)
    fill_nr(db, 'nrseqs')

    db.close()
