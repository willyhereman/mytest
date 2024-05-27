
/* sp_5kdv.com */

batch("sp_setup.max");
batch("sp_5kdv.dat");
batch("sp_elist.max");
writefile("sp_5kdv.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_5kdv.com */
