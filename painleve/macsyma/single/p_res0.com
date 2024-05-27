
batch("np_setup.max");
batch("p_res0.dat");
writefile("p_res0.out");
batch("np_exec.max");
closefile();
/* quit(); */

