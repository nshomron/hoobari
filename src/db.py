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
            res = self.con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='variants'")

            # Drop existing database if needed
            if res.fetchone():
                self.con.execute('DROP TABLE IF EXISTS variants')
                self.con.execute('DROP TABLE IF EXISTS qnames')
                self.con.execute('DROP VIEW IF EXISTS shared_lengths')
                self.con.execute('DROP VIEW IF EXISTS fetal_lengths')

            if not res.fetchone():

                #TODO: Are properties like chromosome constant for each qname
                self.con.execute(   '''CREATE TABLE qnames(
                                    qname varchar(50) PRIMARY KEY,
                                    length int(10) DEFAULT NULL,
                                    is_fetal tinyint(1) DEFAULT NULL
                                    )''')
                self.con.execute('CREATE INDEX idx_length_fetal ON qnames (length, is_fetal)')
                self.con.execute(   '''CREATE TABLE variants(
                                    qname varchar(50) NOT NULL,
                                    chromosome char(2) DEFAULT NULL,
                                    pos int(10) NOT NULL,
                                    genotype varchar(50) DEFAULT NULL,
                                    var_type tinyint(1) DEFAULT NULL,
                                    for_ff tinyint(1) DEFAULT NULL,
                                    FOREIGN KEY(qname) REFERENCES qnames(qname)
                                    )''')
                self.con.execute('CREATE INDEX idx_chrom_pos ON variants (chromosome, pos)')


    # Insert variants to table
    def insertVariant(self, chromosome, position, info_list):

        # Insert qname if it's not in already
        for line in info_list:

            # If the qname is marked as fetal, update accordingly
            if line[3] == 1:
                self.con.execute('''
                INSERT OR REPLACE INTO qnames (length, qname, is_fetal)
                VALUES(:len, ":qname", :fetal)
                ''',{"len":int(line[1]), "qname":line[2], "fetal":int(line[3])})
            else:
                self.con.execute('''
                INSERT OR IGNORE INTO qnames (length, qname, is_fetal)
                VALUES(:len, ":qname", :fetal)
                ''',{"len":int(line[1]), "qname":line[2], "fetal":int(line[3])})



        query = '''
            INSERT INTO `variants`
            (qname,
            chromosome,
            pos,
            genotype,
            var_type,
            for_ff)
            VALUES
            '''


        for line in info_list:
            query += '("{0}","{1}",{2},"{3}",{4},{5}),'.format(line[2],
                    chromosome, position, line[0],  line[4], line[5])

        query = query[:-1]
        self.con.execute(query)
        self.con.commit()

    # Create views for both shared and fetal lengths
    def lengthDists(self):
        self.con.execute('''
        CREATE VIEW fetal_lengths AS SELECT length, qname FROM qnames WHERE is_fetal=1
        ''')
        self.con.execute('''
        CREATE VIEW shared_lengths AS SELECT length, qname FROM qnames WHERE is_fetal=0
        ''')


    # Get fetal lengths
    def getFetalLengths(self):
        return pd.read_sql_query("""
                                    SELECT length, COUNT(length)
                                    FROM fetal_lengths
                                    WHERE chromosome not in ('X','Y')
                                    GROUP BY length
                                """)

    # Get shared lengths
    def getSharedLengths(self):
        return pd.read_sql_query("""
                                    SELECT length, COUNT(length)
                                    FROM shared_lengths
                                    WHERE chromosome not in ('X','Y')
                                    GROUP BY length
                                """)

    # Gets fetal and shared qnames
    def getFetalSharedQnames(self):
        return set(this.con.execute("SELECT DISTINCT(qname) FROM fetal_lengths")), set(this.con.execute("SELECT DISTINCT(qname) FROM shared_lengths"))

    # Gets all variants in specified chromosomal position
    def getPositionVariants(self, chromosome, position):
        return this.con.execute("""
                            SELECT v.genotype, q.length, q.is_fetal
                            FROM variants v, qnames q
                            WHERE v.chromosome=":chr" AND v.position=:pos
                """,{"chr":chromosome, "pos":position})
