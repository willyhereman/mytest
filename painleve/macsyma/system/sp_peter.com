
/* sp_peter.com */

batch("sp_setup.max");
batch("sp_peter.dat");
batch("sp_elist.max");
writefile("sp_peter.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_peter.com */
