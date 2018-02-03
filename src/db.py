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
                self.con.execute('DROP VIEW IF EXISTS shared_length_counts')
                self.con.execute('DROP VIEW IF EXISTS fetal_length_counts')

            if not res.fetchone():

                #TODO: Are properties like chromosome constant for each qname
                self.con.execute(   '''CREATE TABLE qnames(
                                    qname varchar(50) PRIMARY KEY,
                                    length int(10) DEFAULT NULL,
                                    for_ff tinyint(1) DEFAULT NULL,
                                    is_fetal tinyint(1) DEFAULT NULL
                                    )''')
                self.con.execute('CREATE INDEX idx_length_fetal ON qnames (length, for_ff)')
                self.con.execute(   '''CREATE TABLE variants(
                                    qname varchar(50) NOT NULL,
                                    chromosome char(2) DEFAULT NULL,
                                    pos int(10) NOT NULL,
                                    genotype varchar(50) DEFAULT NULL,
                                    var_type tinyint(1) DEFAULT NULL,
                                    FOREIGN KEY(qname) REFERENCES qnames(qname)
                                    )''')
                self.con.execute('CREATE INDEX idx_chrom_pos ON variants (chromosome, pos)')


    # Insert variants to table
    def insertVariant(self, chromosome, position, info_list):

        # Insert qname if it's not in already
        if chromosome not in ('M','X','Y'):
            for line in info_list:

                # If the qname is marked as for_ff, update accordingly
                if line[5] in (1,2):
                    self.con.execute('''
                    INSERT OR REPLACE INTO qnames (length, qname, for_ff, is_fetal)
                    VALUES(:len, :qname, :for_ff, :is_fetal)
                    ''',{"len":int(line[1]), "qname":str(line[2]), "for_ff":int(line[5]), "is_fetal":int(line[3])})
                else:
                    self.con.execute('''
                    INSERT OR IGNORE INTO qnames (length, qname, for_ff, is_fetal)
                    VALUES(:len, :qname, :for_ff, :is_fetal)
                    ''',{"len":int(line[1]), "qname":str(line[2]), "for_ff":int(line[5]), "is_fetal":int(line[3])})



        query = '''
            INSERT INTO `variants`
            (qname,
            chromosome,
            pos,
            genotype,
            var_type)
            VALUES
            '''


        for line in info_list:
            query += '("{0}","{1}",{2},"{3}",{4}),'.format(line[2],
                    chromosome, position, line[0], line[4])

        query = query[:-1]
        self.con.execute(query)
        self.con.commit()

    # Create views for both shared and fetal lengths
    def lengthDists(self):
        self.con.execute('''
        CREATE VIEW fetal_lengths AS SELECT length, qname FROM qnames WHERE for_ff=1
        ''')
        self.con.execute("""
                        CREATE VIEW fetal_length_counts
                        AS SELECT length, COUNT(length)
                        FROM fetal_lengths
                        GROUP BY length
                        """)
        self.con.execute('''
        CREATE VIEW shared_lengths AS SELECT length, qname FROM qnames WHERE for_ff=2
        ''')
        self.con.execute("""
                        CREATE VIEW shared_length_counts
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
    def getPositionVariants(self, chromosome, position):
        return self.con.execute("""
                            SELECT v.genotype, q.length, q.for_ff
                            FROM variants v, qnames q
                            WHERE v.chromosome=":chr" AND v.pos=:pos
                """,{"chr":chromosome, "pos":position})
