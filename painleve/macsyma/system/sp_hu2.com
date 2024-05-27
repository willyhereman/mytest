
/* sp_hu2.com */

batch("sp_setup.max");
batch("sp_hu2.dat");
batch("sp_elist.max");
writefile("sp_hu2.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_hu2.com */
