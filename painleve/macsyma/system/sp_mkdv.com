
/* sp_mkdv.com */

batch("sp_setup.max");
batch("sp_mkdv.dat");
batch("sp_elist.max");
writefile("sp_mkdv.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_mkdv.com */
