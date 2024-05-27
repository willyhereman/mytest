
batch("np_setup.max");
batch("p_haj.dat");
writefile("p_haj.out");
batch("np_exec.max");
closefile();
/* quit(); */

