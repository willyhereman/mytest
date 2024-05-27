
batch("np_setup.max");
batch("p_5pkdv.dat");
writefile("p_5pkdv.out");
batch("np_exec.max");
closefile();
/* quit(); */

