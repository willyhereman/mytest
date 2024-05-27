
/* sp_wave.com */

batch("sp_setup.max");
batch("sp_wave.dat");
batch("sp_elist.max");
writefile("sp_wave.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_wave.com */
