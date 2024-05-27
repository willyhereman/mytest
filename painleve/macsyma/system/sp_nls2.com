
/* sp_nls2.com */

batch("sp_setup.max");
batch("sp_nls2.dat");
batch("sp_elist.max");  
writefile("sp_nls2.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_nls2.com */
