import sqlite3
import sys

# Get database names
dbs=[]
with open(sys.argv[1]) as dbsfile:
    for line in dbsfile:
        dbs.append(line)
db1=dbs[0]

# Merge dbs
for db2 in dbs[1:]:
    MergeDBs(db1,db2)

# Merge the two sqlite DBs with no duplications
def MergeDBs(db1, db2, table='variants'):
    con = sqlite3.connect(db1)

    con.execute("""attach '{0}' as toMerge;
                BEGIN;
                insert or ignore into {1} select * from toMerge.{1} where not exists (select * from {1} where chromosome=toMerge.chromosome AND pos=toMerge.pos);
                COMMIT;
                detach toMerge;""".format(db2,table))
