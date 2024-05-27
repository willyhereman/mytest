
batch("np_setup.max");
batch("p_mkdv.dat");
writefile("p_mkdv.out");
batch("np_exec.max");
closefile();
/* quit(); */

