
batch("np_setup.max");
batch("p_kin.dat");
writefile("p_kin.out");
batch("np_exec.max");
closefile();
/* quit(); */

