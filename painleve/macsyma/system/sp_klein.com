
/* sp_klein.com */

batch("sp_setup.max");
batch("sp_klein.dat");
batch("sp_elist.max");
writefile("sp_klein.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_klein.com */
