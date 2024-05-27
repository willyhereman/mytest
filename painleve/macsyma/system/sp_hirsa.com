
/* sp_hirsa.com */

batch("sp_setup.max");
batch("sp_hirsa.dat");
batch("sp_elist.max");  
writefile("sp_hirsa.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_hirsa.com */
