
batch("np_setup.max")$
batch("p_ks.dat");
writefile("p_ks.out");
batch("np_exec.max")$
closefile()$
/* quit(); */

