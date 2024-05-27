
batch("np_setup.max");
batch("p_tony.dat");
writefile("p_tony.out");
batch("np_exec.max");
closefile();
/* quit(); */

