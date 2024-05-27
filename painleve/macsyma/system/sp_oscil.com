
/* sp_oscil.com */

batch("sp_setup.max");
batch("sp_oscil.dat");
batch("sp_elist.max");
writefile("sp_oscil.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_oscil.com */
