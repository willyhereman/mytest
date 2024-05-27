
batch("np_setup.max");
batch("p_cylkdv.dat");
writefile("p_cylkdv.out");
batch("np_exec.max");
closefile();

