
/* sp_heno1.com */

batch("sp_setup.max");
batch("sp_heno1.dat");
batch("sp_elist.max");
writefile("sp_heno1.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_heno1.com */
