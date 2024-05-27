
/* sp_7kdv.com */

batch("sp_setup.max");
batch("sp_7kdv.dat"); 
batch("sp_elist.max");
writefile("sp_7kdv.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_7kdv.com */
