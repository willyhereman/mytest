
batch("np_setup.max");
batch("p_wil.dat");
writefile("p_wil.out");
batch("np_exec.max");
closefile();
/* quit(); */

