
/* sp_stars.com */

batch("sp_setup.max");
batch("sp_stars.dat");
batch("sp_elist.max");
writefile("sp_stars.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_stars.com */
