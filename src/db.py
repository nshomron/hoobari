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
            
            self.con.execute('CREATE INDEX idx_chrom_pos ON variants (chromosome, pos);')

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
<<<<<<< 9fb2af6a7dd3367ee9066cb663ffac6ed6e9267e
        return pd.read_sql_query("select distinct(`length`) as len, count(*) as `count` from variants where for_ff=2 and chromosome not in ('X', 'Y') group by len", self.con)
=======
        return pd.read_sql_query("select * from shared_lengths", self.con)

    # Create length distribution table
    def createDistTable(self):
        self.con.execute('''
        create table fetal_lengths(
            `length` int NOT NULL,
            `count`  int NOT NULL DEFAULT '0'
            PRIMARY KEY (`length`),
            UNIQUE KEY `length` (`length`)
        )
        ''')

        self.con.execute('''
        create table shared_lengths(
            `length` int(5) NOT NULL,
            `count`  int(5) NOT NULL DEFAULT '0'
            PRIMARY KEY (`length`),
            UNIQUE KEY `length` (`length`)
        )
        ''')
>>>>>>> Removed unsupported unsigned type

    def update_is_fetal(self):
        self.con.execute('UPDATE variants SET is_fetal=1 WHERE qname=(SELECT qname FROM variants WHERE is_fetal=1)')
