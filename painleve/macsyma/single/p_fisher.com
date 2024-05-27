
/* Command file for Fisher equation in ODE form */

/* ************************************************************************ */
/*                          Batch file P_FISHER.COM                         */
/* ************************************************************************ */

batch("np_setup.max")$
batch("p_fisher.dat")$
writefile("p_fisher.out");
batch("np_exec.max")$
closefile();
quit();

/* ************************** END of P_FISHER.COM ************************* */

