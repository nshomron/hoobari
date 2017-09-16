create table if not exists fetal_lengths(
    `length` int NOT NULL,
    `count`  int NOT NULL DEFAULT '0',
    PRIMARY KEY (`length`)
<<<<<<< HEAD
);
=======
)
>>>>>>> 9b2c6c108b08a44af3c61a01a42a49d26c498d7d

create table if not exists shared_lengths(
    `length` int NOT NULL,
    `count`  int NOT NULL DEFAULT '0',
    PRIMARY KEY (`length`)
<<<<<<< HEAD
);
=======
)
>>>>>>> 9b2c6c108b08a44af3c61a01a42a49d26c498d7d

insert into fetal_lengths 
select `length`, count(*) as `count` 
from (select min(`length`) as `length` 
      from variants where for_ff=1 and chromosome not in ('X', 'Y') 
      group by `qname`) as qunique 
<<<<<<< HEAD
group by `length`;
=======
group by `length`
>>>>>>> 9b2c6c108b08a44af3c61a01a42a49d26c498d7d

insert into shared_lengths 
select `length`, count(*) as `count` 
from (select min(`length`) as `length` 
      from variants where for_ff=2 and chromosome not in ('X', 'Y') 
      group by `qname`) as qunique 
<<<<<<< HEAD
group by `length`;
=======
group by `length`
>>>>>>> 9b2c6c108b08a44af3c61a01a42a49d26c498d7d
