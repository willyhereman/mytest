
/* sp_redif.com */

batch("sp_setup.max");
batch("sp_redif.dat");
batch("sp_elist.max");
writefile("sp_redif.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_redif.com */
