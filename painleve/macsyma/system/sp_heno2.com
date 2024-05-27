
/* sp_heno2.com */

batch("sp_setup.max");
batch("sp_heno2.dat");
batch("sp_elist.max");
writefile("sp_heno2.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_heno2.com */
