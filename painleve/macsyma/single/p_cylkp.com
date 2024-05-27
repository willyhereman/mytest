
/* This is command file for the cylindrical KP equation */
/* in book Ablowitz and Clarkson on p 64 */

batchload("np_setup.max");
writefile("p_cylkp.out");
batch("p_cylkp.dat");
batch("np_exec.max");
closefile();
/* quit(); */

