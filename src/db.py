import sqlite3
import csv
from sys import stderr

class Variants(object):
    def __init__(self, dropdb=False,dbname='hoobari.db'):
        # Connect to DB
        self.con=sqlite3.connect(dbname,isolation_level=None)

        # Drop existing database if needed
        if dropdb:
            self.con.execute('DROP TABLE `variants`;')

        # Check table's existance
        res=self.con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='variants';")
        if not res.fetchone():
            self.con.execute('CREATE TABLE `variants` (`pos` int(10) NOT NULL,`chromosome` char(2) DEFAULT NULL,`genotype` char(1) DEFAULT NULL,`qname` varchar(50) DEFAULT NULL)')

    # Insert variants to table
    def insertVariants(self,genotypes, chromosome, position):
        query='''
            INSERT INTO `variants`
            (`pos`,
            `chromosome`,
            `genotype`,
            `qname`,
            `length`)
            VALUES
            '''
        for genotype in genotypes:
            query+='({0},"{1}","{2}","{3}",{4}),'.format(position,
                    chromosome,genotype[0],genotype[2],genotype[1])
        query=query[:-1] + ';'
        self.con.execute(query)

    # Insert variant to table
    def insertVariant(length,pos,chromosome,genotype,qname):
        self.con.execute('''
            INSERT INTO `variants`
            (`pos`,
            `chromosome`,
            `genotype`,
            `qname`,
            `length`)
            VALUES
            ({0},
            {1},
            {2},
            {3},
            {4});
        '''.format(pos,chromosomes,genotype,qname,length))
