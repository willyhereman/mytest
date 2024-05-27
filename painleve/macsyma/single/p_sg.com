
/* Command file for the Painleve test of the sine-Gordon equation */

/* ************************************************************************* */
/*                           Batch file P_SG.COM                             */
/* ************************************************************************* */

batch("np_setup.max")$
batch("p_sg.dat")$
writefile("p_sg.out");
batch("np_exec.max")$
closefile();
quit();

/* **************************** END of P_SG.COM *************************** */

