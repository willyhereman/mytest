
batch("np_setup.max");
batch("p_gil.dat");
writefile("p_gil.out");
batch("np_exec.max");
closefile();
/* quit(); */

