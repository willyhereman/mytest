
batch("np_setup.max");
batch("p_varco.dat");
writefile("p_varco.out");
batch("np_exec.max");
closefile();
/* quit(); */

