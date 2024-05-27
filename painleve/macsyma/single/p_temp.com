
/* Command file for 1 dem dipole equation in ODE form */

/* ************************************************************************ */
/*                          Batch file P_DIPOLE.COM                         */
/* ************************************************************************ */

batch("np_setup.max")$
batch("p_dipole.dat")$
writefile("p_dipole.out");
batch("np_exec.max")$
closefile();
quit();

/* ************************** END of P_DIPOLE.COM ************************* */

