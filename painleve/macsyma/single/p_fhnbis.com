
batch("np_setup.max");
batch("p_fhnbis.dat");
writefile("p_fhnbis.out");
batch("np_exec.max");
closefile();
/* quit(); */

