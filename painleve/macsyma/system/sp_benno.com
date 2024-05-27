
/* sp_benno.com */

batch("sp_setup.max");
batch("sp_benno.dat");
batch("sp_elist.max");
writefile("sp_benno.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_benno.com */
