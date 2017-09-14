update variants set for_ff=1 where qname in (select distinct(qname) from variants where for_ff=1);
update variants set for_ff=2 where qname in (select distinct(qname) from variants where for_ff=2);
