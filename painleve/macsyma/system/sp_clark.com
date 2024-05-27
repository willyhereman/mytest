
/* sp_clark.com */

batch("sp_setup.max");
batch("sp_clark.dat");
batch("sp_elist.max");
writefile("sp_clark.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_clark.com */
