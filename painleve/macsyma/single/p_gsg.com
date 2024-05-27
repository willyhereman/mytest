
batch("np_setup.max");
batch("p_gsg.dat");
writefile("p_gsg.out");
batch("np_exec.max");
closefile();
/* quit(); */

