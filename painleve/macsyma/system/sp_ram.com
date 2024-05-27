
/* sp_ram.com */

batch("sp_setup.max");
batch("sp_ram.dat"); 
batch("sp_elist.max");
writefile("sp_ram.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_ram.com */
