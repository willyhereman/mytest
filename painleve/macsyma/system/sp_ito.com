
/* sp_ito.com */

batch("sp_setup.max");
batch("sp_ito.dat");
batch("sp_elist.max");
writefile("sp_ito.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_ito.com */
