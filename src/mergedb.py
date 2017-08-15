import sqlite3

# Merge the two sqlite DBs with no duplications
def MergeDBs(db1, db2, table):
    con = sqlite3.connect(db1)

    con.execute("""attach '{0}' as toMerge;
                BEGIN;
                insert or ignore into {1} select * from toMerge.{1};
                COMMIT;
                detach toMerge;""".format(db2,table))
