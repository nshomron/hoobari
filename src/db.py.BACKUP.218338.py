import sqlite3
import os
import pandas as pd

class Variants(object):
    def __init__(self, dbpath = './hoobari.db', probe=True):

        # create directory
        os.makedirs(os.path.dirname(dbpath), exist_ok=True)

        # Connect to DB
        self.con = sqlite3.connect(dbpath, isolation_level = None)

        # Check table's existance
<<<<<<< HEAD
        res = self.con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='variants';")

        # Drop existing database if needed
        if res.fetchone():
            self.con.execute('DROP TABLE `variants`;')

        if not res.fetchone():
            self.con.execute(   '''CREATE TABLE IF NOT EXISTS `variants`(
                                `chromosome` char(2) DEFAULT NULL,
                                `pos` int(10) NOT NULL,
                                `genotype` varchar(50) DEFAULT NULL,
                                `length` int(10) DEFAULT NULL,
                                `qname` varchar(50) DEFAULT NULL,
                                `is_fetal` tinyint(1) DEFAULT NULL,
                                `var_type` tinyint(1) DEFAULT NULL,
                                `for_ff` tinyint(1) DEFAULT NULL)''')
            
            self.con.execute('CREATE INDEX idx_chrom_pos ON variants (chromosome, pos);')
=======
        if probe:
            res = self.con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='variants';")

            # Drop existing database if needed
            if res.fetchone():
                self.con.execute('DROP TABLE `variants`;')

            if not res.fetchone():
                self.con.execute(   '''CREATE TABLE IF NOT EXISTS `variants`(
                                    `chromosome` char(2) DEFAULT NULL,
                                    `pos` int(10) NOT NULL,
                                    `genotype` varchar(50) DEFAULT NULL,
                                    `length` int(10) DEFAULT NULL,
                                    `qname` varchar(50) DEFAULT NULL,
                                    `is_fetal` tinyint(1) DEFAULT NULL,
                                    `var_type` tinyint(1) DEFAULT NULL,
                                    `for_ff` tinyint(1) DEFAULT NULL)''')
>>>>>>> Better length distribution direct from database

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
        self.con.execute(query) ### TODO: check if execute many is better
        #self.con.commit()

    def fetalLengthDist(self):
        return pd.read_sql_query("select distinct(`length`) as len, count(*) as `count` from variants where for_ff=1 and chromosome not in ('X', 'Y') group by len", self.con)

    def sharedLengthDist(self):
        return pd.read_sql_query("select distinct(`length`) as len, count(*) as `count` from variants where for_ff=2 and chromosome not in ('X', 'Y') group by len", self.con)
