
batch("np_setup.max");
batch("p_ram2.dat");
writefile("p_ram2.out");
batch("np_exec.max");
closefile();
/* quit(); */

