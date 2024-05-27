
/* sp_nls1.com */

batch("sp_setup.max");
batch("sp_nls1.dat");
batch("sp_elist.max");  
writefile("sp_nls1.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_nls1.com */
