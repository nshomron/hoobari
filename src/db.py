import sqlite3


class Variants(Object):
    def __init__(self,dbname='hoobari.db'):
        # Connect to DB
        con=sqlite3.connect(dbname,isolation_level=None)

        # Check table's existance
        res=con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='variants';")
        if not res.fetchone():
            con.execute('CREATE TABLE `variants` (`pos` int(10) NOT NULL,`chromosome` char(2) DEFAULT NULL,`genotype` char(1) DEFAULT NULL,`qname` varchar(50) DEFAULT NULL)')

    # Insert variants to table
    def insertVariant(genotypes, chromosome, position):
        query='''
            INSERT INTO `variants`
            (`pos`,
            `chromosome`,
            `genotype`,
            `qname`)
            VALUES
            (
            '''
        for genotype in genotypes:
            query+='({0},{1},{2},{3}),'.format(pos,chromosomes,genotype,qname)
        con.execute(query[:-1] + ');')

    # Insert variant to table
    def insertVariant(pos,chromosome,genotype,qname)
        con.execute('''
            INSERT INTO `variants`
            (`pos`,
            `chromosome`,
            `genotype`,
            `qname`)
            VALUES
            ({0},
            {1},
            {2},
            {3});
        '''.format(pos,chromosomes,genotype,qname))
