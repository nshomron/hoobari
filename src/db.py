import sqlite3
import csv
from sys import stderr

class Variants(object):
    def __init__(self, dropdb=False, dbpath='.'):
        # Connect to DB
        self.con = sqlite3.connect(dbpath + '/hoobari.db', isolation_level = None, timeout=1000)

        # Drop existing database if needed
        if dropdb:
            self.con.execute('DROP TABLE `variants`;')

        # Check table's existance
        res = self.con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='variants';")
        
        if not res.fetchone():
            try:
                self.con.execute('CREATE TABLE `variants` (`chromosome` char(2) DEFAULT NULL,`pos` int(10) NOT NULL,`genotype` char(20) DEFAULT NULL,`length` int(10) DEFAULT NULL,`qname` varchar(50) DEFAULT NULL)')
            except sqlite3.OperationalError:
                pass

    # Insert variants to table
    def insertVariant(self, chromosome, position, info_list):
        
        query = '''
            INSERT INTO `variants`
            (`chromosome`,
            `pos`,
            `genotype`,
            `length`,
            `qname`)
            VALUES
            '''
        
        for line in info_list:
            query += '("{0}",{1},"{2}",{3},"{4}"),'.format(chromosome,
                        position, line[0], line[1], line[2])
        
        query = query[:-1] + ';'
        
        self.con.execute(query)
        self.con.commit()