
batch("np_setup.max");
batch("p_bous.dat");
writefile("p_bous.out");
batch("np_exec.max");
closefile();
/* quit(); */

