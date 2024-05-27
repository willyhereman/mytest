
/* sp_sine.com */

batch("sp_setup.max");
batch("sp_sine.dat");
batch("sp_elist.max");
writefile("sp_sine.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_sine.com */
