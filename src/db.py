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
            # Drop existing database if needed
            res_list = [i[0:2] for i in self.con.execute("SELECT type, name FROM sqlite_master WHERE type in ('table','view')").fetchall()]
            for res in res_list:
                self.con.execute('DROP ' + res[0] + ' IF EXISTS ' + res[1])

            #TODO: Are properties like chromosome constant for each qname
            self.con.execute(   '''CREATE TABLE qnames(
                                qname varchar(50) PRIMARY KEY,
                                is_fetal tinyint(1) DEFAULT NULL
                                )''')
            # self.con.execute('CREATE INDEX idx_length_fetal ON qnames (length, for_ff)')
            self.con.execute(   '''CREATE TABLE variants(
                                qname varchar(50) NOT NULL,
                                chromosome char(2) DEFAULT NULL,
                                pos int(10) NOT NULL,
                                genotype varchar(50) DEFAULT NULL,
                                var_type tinyint(1) DEFAULT NULL,
                                length int(10) DEFAULT NULL,
                                for_ff tinyint(1) DEFAULT NULL,
                                is_fetal tinyint(1) DEFAULT NULL,
                                FOREIGN KEY(qname) REFERENCES qnames(qname)
                                )''')
            self.con.execute('CREATE INDEX idx_chrom_pos ON variants (chromosome, pos)')
            self.con.execute(   '''CREATE TABLE samples(
                                mother char(20) DEFAULT NULL,
                                father char(20) DEFAULT NULL)''')

    # Insert variants to table
    def insertVariant(self, chromosome, position, info_list):

        query = '''
            INSERT INTO `variants`
            (qname,
            chromosome,
            pos,
            genotype,
            var_type,
            length,
            for_ff,
            is_fetal)
            VALUES
            '''


        for line in info_list:
            query += '("{0}","{1}",{2},"{3}",{4},{5},{6},{7}),'.format(line[2],
                    chromosome, position, line[0], line[4], line[1], line[5], line[3])

        query = query[:-1]
        self.con.execute(query)
        self.con.commit()

    # Create qnames table for origin model
    def createQnamesTable(self):
        self.con.execute('''INSERT INTO qnames 
                                    SELECT qname, max(is_fetal) 
                                    from variants 
                                    group by qname''')

    # Create views for both shared and fetal lengths
    def lengthDists(self):
        self.con.execute('''
        CREATE VIEW fetal_lengths AS SELECT length, qname FROM variants WHERE for_ff=1 GROUP BY qname
        ''')
        self.con.execute("""
                        CREATE TABLE fetal_length_counts
                        AS SELECT length, COUNT(length)
                        FROM fetal_lengths
                        GROUP BY length
                        """)
        self.con.execute('''
        CREATE VIEW shared_lengths AS SELECT length, qname FROM variants WHERE for_ff=2 GROUP BY qname
        ''')
        self.con.execute("""
                        CREATE TABLE shared_length_counts
                        AS SELECT length, COUNT(length)
                        FROM shared_lengths
                        GROUP BY length
                        """)

    # Get fetal lengths
    def getFetalLengths(self):
        return pd.read_sql_query("select * from fetal_length_counts", self.con).set_index('length')

    # Get shared lengths
    def getSharedLengths(self):
        return pd.read_sql_query("select * from shared_length_counts", self.con).set_index('length')

    # Gets fetal and shared qnames
    def getFetalSharedQnames(self):
        return set(self.con.execute("SELECT DISTINCT(qname) FROM fetal_lengths")), set(self.con.execute("SELECT DISTINCT(qname) FROM shared_lengths"))

    # Gets all variants in specified chromosomal position
    def getPositionVariants(self, chromosome, position, model):
        if model == 'origin':
            return self.con.execute("""
                                SELECT v.genotype, v.length, q.is_fetal
                                FROM variants v, qnames q
                                WHERE v.chromosome=:chr AND v.pos=:pos AND q.qname=v.qname
                    """,{"chr":chromosome, "pos":position})
        else:
            return self.con.execute("""
                                SELECT genotype, length, is_fetal
                                FROM variants
                                WHERE chromosome=:chr AND pos=:pos
                    """,{"chr":chromosome, "pos":position})

    def create_samples_table(self, parental_samples):
        mother, father = parental_samples
        self.con.execute('INSERT INTO `samples` (mother, father) VALUES (:m, :f)', {'m': str(mother), 'f': str(father)})

    def get_samples(self):
        return self.con.execute('SELECT * FROM samples').fetchall()[0]
