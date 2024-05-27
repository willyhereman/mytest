
/* sp_kaup.com */

batch("sp_setup.max");
batch("sp_kaup.dat");
batch("sp_elist.max");
writefile("sp_kaup.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_kaup.com */
