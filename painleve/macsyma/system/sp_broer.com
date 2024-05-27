
/* sp_broer.com */

batch("sp_setup.max");
batch("sp_broer.dat");
batch("sp_elist.max");
writefile("sp_broer.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_broer.com */
