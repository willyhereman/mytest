
/* Command file for 1 dem dipol2 equation in ODE form */

/* ************************************************************************ */
/*                          Batch file P_dipol2.COM                         */
/* ************************************************************************ */

batch("np_setup.max")$
batch("p_dipol2.dat")$
writefile("p_dipol2.out");
batch("np_exec.max")$
closefile();
quit();

/* ************************** END of P_dipol2.COM ************************* */

