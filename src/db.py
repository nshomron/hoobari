import pymysql as mysql
import csv
from sys import stderr
from re import split as resplit

class Variants(object):
    def __init__(self, mysql_info, db_name, dropdb=False):
        # Connect to DB
        self.con = self.ConnectDB(mysql_info, db_name)
        self.cur = self.con.cursor()

        # Create database if there isn't one already
        self.cur.execute('CREATE DATABASE IF NOT EXISTS {0}'.format(db_name))
        self.con.select_db(db_name)

        # Drop existing table if needed
        if dropdb:
            self.cur.execute('DROP TABLE `variants`')

        # Create table if there isn't one already
        res = self.cur.execute("""
            CREATE TABLE IF NOT EXISTS `variants` (
            `chromosome` char(2) DEFAULT NULL,
            `pos` int(10) unsigned NOT NULL,
            `genotype` varchar(50) DEFAULT NULL,
            `length` int(10) DEFAULT NULL,
            `qname` varchar(50) DEFAULT NULL,
            `is_fetal` tinyint(1) DEFAULT NULL,
            `var_type` tinyint(1) DEFAULT NULL,
            `for_ff` tinyint(1) DEFAULT NULL
            )
        """)

        # Commit changes
        self.con.commit()

    def ConnectDB(self, mysql_info, db_name):
    
        first_split = mysql_info.split(':')
        port = int(first_split[1])
        user, host, socket = resplit(r'@|/', first_split[0], 2)
        socket = '/' + socket
        
        connection = mysql.connect(host = host, user = user, port = port, database = db_name, unix_socket = socket)

        return connection

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
            query += '("{0}",{1},"{2}",{3},"{4}",{5},{6},{7}),'.format(chromosome, position, line[0], line[1], line[2], line[3], line[4], line[5])

        query = query[:-1] + ';'

        self.cur.execute(query)
        self.con.commit()




# import sqlite3
# from sys import stderr
# import os

# class Variants(object):
#     def __init__(self, dropdb=False, dbpath='.'):
#         # Connect to DB
#         self.con = sqlite3.connect(dbpath, isolation_level = None)

#         # Check table's existance
#         res = self.con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='variants';")

#         # Drop existing database if needed
#         if dropdb and res.fetchone():
#             self.con.execute('DROP TABLE `variants`;')

#         if not res.fetchone():
#             self.con.execute(   'CREATE TABLE `variants` (' +
#                                 '`chromosome` char(2) DEFAULT NULL,' + 
#                                 '`pos` int(10) NOT NULL,' +
#                                 '`genotype` char(20) DEFAULT NULL,' + 
#                                 '`length` int(10) DEFAULT NULL,' +
#                                 '`qname` varchar(50) DEFAULT NULL,' + 
#                                 '`is_fetal` int(1) DEFAULT NULL,' +
#                                 '`var_type` numeric(6) DEFAULT NULL,' +
#                                 '`for_ff` int(1) DEFAULT NULL)')

#     # Insert variants to table
#     def insertVariant(self, chromosome, position, info_list):
        
#         query = '''
#             INSERT INTO `variants`
#             (`chromosome`,
#             `pos`,
#             `genotype`,
#             `length`,
#             `qname`,
#             `is_fetal`,
#             `var_type`,
#             `for_ff`)
#             VALUES
#             '''
        
#         for line in info_list:
#             query += '("{0}",{1},"{2}",{3},"{4}",{5},{6},{7}),'.format(chromosome,
#                         position, line[0], line[1], line[2], line[3], line[4], line[5])
        
#         query = query[:-1] + ';'
        
#         self.con.execute(query)