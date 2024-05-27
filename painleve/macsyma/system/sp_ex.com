
/* sp_ex.com */

batch("sp_setup.max");
batch("sp_ex.dat");
batch("sp_elist.max");
writefile("sp_ex.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_ex.com */
