
/* Command file for a simple ODE */

/* ************************************************************************ */
/*                          Batch file P_SIMODE.COM                         */
/* ************************************************************************ */

batch("np_setup.max")$
batch("p_simode.dat")$
writefile("p_simode.out");
batch("np_exec.max")$
closefile();
quit();

/* ************************** END of P_SIMODE.COM ************************* */

