
/* sp_univ.com */

batch("sp_setup.max");
batch("sp_univ.dat");
batch("sp_elist.max"); 
writefile("sp_univ.out");
batch("sp_exec.max");
closefile();
/* quit(); */ 

/* sp_univ.com */
