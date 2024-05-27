
/* sp_bur.com */

batch("sp_setup.max");
batch("sp_bur.dat");
batch("sp_elist.max");
writefile("sp_bur.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_bur.com */
