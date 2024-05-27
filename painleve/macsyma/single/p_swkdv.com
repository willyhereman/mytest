
batch("np_setup.max");
batch("p_swkdv.dat");
writefile("p_swkdv.out");
batch("np_exec.max");
closefile();
/* quit(); */

