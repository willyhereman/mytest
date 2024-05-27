
/* Command file for the Painleve test of the alternate sine-Gordon equation */

/* ************************************************************************* */
/*                           Batch file P_ASG.COM                            */
/* ************************************************************************* */

batch("np_setup.max")$
batch("p_asg.dat")$
writefile("p_asg.out");
batch("np_exec.max")$
closefile();
quit();

/* **************************** END of P_ASG.COM ************************** */

