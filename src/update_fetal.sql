update variants set is_fetal=1 where qname in (select qname from variants where is_fetal=1);
