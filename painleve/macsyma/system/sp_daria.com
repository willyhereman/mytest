
/* sp_daria.com */

batch("sp_setup.max");
writefile("sp_daria.out");
batch("sp_daria.dat");
batch("sp_elist.max");  
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_daria.com */

