
/* Command file for 1 dem dipol1 equation in ODE form */

/* ************************************************************************ */
/*                          Batch file P_dipol1.COM                         */
/* ************************************************************************ */

batch("np_setup.max")$
batch("p_dipol1.dat")$
writefile("p_dipol1.out");
batch("np_exec.max")$
closefile();
quit();

/* ************************** END of P_dipol1.COM ************************* */

