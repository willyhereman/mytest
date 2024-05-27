
batch("np_setup.max");
batch("p_ram1.dat");
writefile("p_ram1.out");
batch("np_exec.max");
closefile();
/* quit(); */

