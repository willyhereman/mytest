
/* sp_hkdv.com */

batch("sp_setup.max");
batch("sp_hkdv.dat");
batch("sp_elist.dat");
writefile("sp_hkdv.out");
batch("sp_exec.max");
closefile();
/* quit(); */

/* sp_hkdv.com */
