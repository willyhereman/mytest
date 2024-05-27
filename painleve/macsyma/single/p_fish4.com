
/* Command file for a fourth-degree Fisher equation in ODE form */

/* ************************************************************************ */
/*                          Batch file P_FIS4.COM                          */
/* ************************************************************************ */

batch("np_setup.max")$
batch("p_fish4.dat")$
writefile("p_fish4.out");
batch("np_exec.max")$
closefile();
/* quit(); */

/* ************************** END of P_FISH4.COM ************************* */

