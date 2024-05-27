
batch("np_setup.max");
batch("p_kdv.dat");
writefile("p_kdv.out");
batch("np_exec.max");
closefile();
/* quit(); */

