
batch("np_setup.max");
batch("p_burg.dat");
writefile("p_burg.out");
batch("np_exec.max");
closefile();
/* quit(); */

