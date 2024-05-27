
/* sp_akak.com */

batch("sp_setup.max");
batch("sp_akak.dat");
batch("sp_elist.max");
writefile("sp_akak.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_akak.com */
