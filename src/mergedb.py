import sqlite3
import sys
import shutil

# Merge the two sqlite DBs with no duplications
def MergeDBs(db1, db2, table='variants'):
    con = sqlite3.connect(db1)

    con.executescript("""attach '{0}' as toMerge;
                         BEGIN;
                         insert or ignore into {1} select * from toMerge.{1}
                         where not exists (select * from {1} where chromosome=toMerge.chromosome AND pos=toMerge.pos);
                         COMMIT;
                         detach toMerge;""".format(db2,table))

# Get database names
dbs=[]
with open(sys.argv[1]) as dbsfile:
    for line in dbsfile:
        dbs.append(line.strip())

# Copy origin database
db1=sys.argv[2]
shutil.copyfile(dbs[0],db1)

# Merge dbs
for db2 in dbs[1:]:
    MergeDBs(db1,db2)

