
/* sp_car.com */

batch("sp_setup.max");
batch("sp_car.dat");
batch("sp_elist.max");
writefile("sp_car.out");
batch("sp_exec.max");
closefile();
/* quit(); */ 

/* sp_car.com */
