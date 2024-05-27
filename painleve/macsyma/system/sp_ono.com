
/* sp_ono.com */

batch("sp_setup.max");
batch("sp_ono.dat");
batch("sp_elist.max");
writefile("sp_ono.out");
batch("sp_exec.max");
closefile();

/* sp_ono.com */
