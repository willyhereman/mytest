
/* sp_so.com */

batch("sp_setup.max");
batch("sp_so.dat");
batch("sp_elist.max");
writefile("sp_so.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_so.com */
