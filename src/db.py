import sqlite3
import csv
from sys import stderr
import os

class Variants(object):
    def __init__(self, dropdb=False, dbpath='.'):
        # Connect to DB
        self.con = sqlite3.connect(dbpath, isolation_level = None)

        # Check table's existance
        res = self.con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='variants';")

        # Drop existing database if needed
        if dropdb and res.fetchone():
            self.con.execute('DROP TABLE `variants`;')

        if not res.fetchone():
            self.con.execute(   'CREATE TABLE `variants` (' +
                                '`chromosome` char(2) DEFAULT NULL,' + 
                                '`pos` int(10) NOT NULL,' +
                                '`genotype` char(20) DEFAULT NULL,' + 
                                '`length` int(10) DEFAULT NULL,' +
                                '`qname` varchar(50) DEFAULT NULL,' + 
                                '`is_fetal` int(1) DEFAULT NULL,' +
                                '`var_type` numeric(6) DEFAULT NULL,' +
                                '`for_ff` int(1) DEFAULT NULL)')

    # Insert variants to table
    def insertVariant(self, chromosome, position, info_list):
        
        query = '''
            INSERT INTO `variants`
            (`chromosome`,
            `pos`,
            `genotype`,
            `length`,
            `qname`,
            `is_fetal`,
            `var_type`,
            `for_ff`)
            VALUES
            '''
        
        for line in info_list:
            query += '("{0}",{1},"{2}",{3},"{4}",{5},{6},{7}),'.format(chromosome,
                        position, line[0], line[1], line[2], line[3], line[4], line[5])
        
        query = query[:-1] + ';'
        
        self.con.execute(query)