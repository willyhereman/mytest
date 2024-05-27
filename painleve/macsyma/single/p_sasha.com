
/* Command file for the Painleve test of the equation due to Sasha Mikhailov */

/* ************************************************************************* */
/*                           Batch file P_SASHA.COM                          */
/* ************************************************************************* */

batch("np_setup.max")$
batch("p_sasha.dat")$
writefile("p_sasha.out");
batch("np_exec.max")$
closefile();
quit();

/* **************************** END of P_SASHA.COM ************************ */

