
batch("np_setup.max");
batch("p_fhn.dat");
writefile("p_fhn.out");
batch("np_exec.max");
closefile();
/* quit(); */ 

