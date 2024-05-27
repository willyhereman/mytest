
/* sp_bous.com */

batch("sp_setup.max");
batch("sp_bous.dat");
batch("sp_elist.max");
writefile("sp_bous.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_bous.com */
