
/* sp_soko.com */

batch("sp_setup.max");
batch("sp_soko.dat");
batch("sp_elist.max");
writefile("sp_soko.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_soko.com */
