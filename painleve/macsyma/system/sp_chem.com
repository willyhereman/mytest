
/* sp_chem.com */

batch("sp_setup.max");
batch("sp_chem.dat");
batch("sp_elist.max");
writefile("sp_chem.out");
batch("sp_exec.412");
closefile();
/* quit(); */

/* sp_chem.com */
