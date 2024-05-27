
batch("np_setup.max");
writefile("p_kru.out");
batch("p_kru.dat");
batch("np_exec.max");
closefile();

