
/* sp_c3ckdv.com */

batch("sp_setup.max");
batch("sp_3ckdv.dat");
batch("sp_elist.max");
writefile("sp_3ckdv.out");
batch("sp_exec.max");
closefile();
/* quit(); */ 

/* sp_c3ckdv.com */
