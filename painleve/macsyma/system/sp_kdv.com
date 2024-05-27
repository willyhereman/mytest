
/* sp_kdv.com */

batch("sp_setup.max");
batch("sp_kdv.dat");
batch("sp_elist.max");
writefile("sp_kdv.out");
batch("sp_exec.max");
closefile();
/* quit(); */ 

/* sp_kdv.com */
