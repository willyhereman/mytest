
/* sp_hu1.com */

batch("sp_setup.max");
batch("sp_hu1.dat");
batch("sp_elist.max");
writefile("sp_hu1.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_hu1.com */
