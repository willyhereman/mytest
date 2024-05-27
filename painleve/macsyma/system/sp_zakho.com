
/* sp_zakho.com */

batch("sp_setup.max");
batch("sp_zakho.dat");
batch("sp_elist.max");
writefile("sp_zakho.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_zakho.com */
