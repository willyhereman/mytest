
batch("np_setup.max");
batch("p_ram3.dat");
writefile("p_ram3.out");
batch("np_exec.max");
closefile();
/* quit(); */

