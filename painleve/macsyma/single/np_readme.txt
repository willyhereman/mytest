
/* ************************************************************************ */
/*                    Documentation file NP_README.TXT                      */
/* ************************************************************************ */

      /**************************************************************/
      /*                                                            */
      /*           *** M A C S Y M A  P R O G R A M ***             */
      /*                                                            */
      /*   THE PAINLEVE TEST FOR SINGLE ORDINARY AND PARTIAL        */
      /*                 DIFFERENTIAL EQUATIONS                     */
      /*                                                            */
      /*  program name: NP_SING.MAX                                 */
      /*  purpose: carry out the Painleve integrability test for    */
      /*           a single ordinary or partial differential eq.    */
      /*           (with or without the Kruskal simplification)     */
      /*  for systems use NP_SYS.MAX which also handles single eqs. */
      /*  computers: tested on VAX-8650, IBM Risc 6000, and         */
      /*             on IBM Compatible PC                           */
      /*  language: MACSYMA release 417.125,                        */
      /*            compatible with Macsyma 2.0 for PC              */
      /*                                                            */
      /*  authors:  W. Hereman, Unal Goktas and Chris Elmer         */
      /*                                                            */
      /*           Department of Mathematical and Computer Sciences */
      /*           Colorado School of Mines,                        */
      /*           Golden, CO 80401-1887, USA                       */
      /*                                                            */
      /*           Version  2.5.  Updated: August 28, 1995          */
      /*                                                            */
      /**************************************************************/

/* *********************************************************************** */
/*                                                                         */
/*     Batch file NP_SING.MAX with the main program for the Painleve test  */
/*                                                                         */
/*              Initial version of the program was written by              */
/*                    Willy HEREMAN and Eric VAN DEN BULCK                 */
/*                                                                         */
/*               Originally developed at: Mathematics Department           */
/*               University of Wisconsin, Madison, WI 53706 USA            */
/*                                                                         */
/*                   Includes the Kruskal simplification                   */
/*                                                                         */
/* *********************************************************************** */

Copyright by Willy Hereman, Eric Van Den Bulck, and Unal Goktas: 

No part of the program NP_SING.MAX for the PAINLEVE TEST of single ODEs 
and PDEs may be reproduced or sold without written consent of the authors.

MACSYMA copyright and trademark: Macsyma, Inc., Arlington, MA, USA.

We are glad to offer you the possibility to carry out the tedious
calculations for the Painleve test for single ODEs and PDEs by computer.

You will be using this version of MACSYMA if your computer has the possibility
to implement Common Lisp Macsyma 412 or higher, including Macsyma 
417.125 and 2.0 for PC.
 
                --------------------------------------

A discussion of the Painleve Test and some other examples may be found in 
the papers:

File macsnews.tex:

W. Hereman and S. Angenent,
``The Painlev\'e test for nonlinear ordinary and partial 
differential equations", MACSYMA Newsletter, vol. 6, pp. 11-8 (1989).

File macsprog.tex:

W. Hereman and E. Van den Bulck,
``MACSYMA program for the Painlev\'e test of nonlinear ordinary 
and partial differential equations", 
Proceedings of the Workshop on Finite Dimensional
Integrable Nonlinear Dynamical Systems, Eds.: P.G.L. Leach and W.-H. Steeb,
Johannesburg, South Africa, January 11-15, 1988, 
World Scientific, Singapore, pp. 117-129 (1988).

More examples of the use of the program in investigation equations with
variable coefficients can be found in:

File kdv95.tex:

W. Hereman and W. Zhuang,
``Symbolic software for soliton theory",
Proceedings of Conference KdV '95, April 1995, Amsterdam, The Netherlands, 
Eds.: M. Hazewinkel, H. W. Capel and E. M. de Jager, 
Kluwer Academic Publishers, pp. 361-378 (1995).
Also: Acta Applicandae Mathematicae, vol 39, pp. 361-378 (1995)

File pdeimacs.tex:

W. Hereman,
``Symbolic software for the study of nonlinear partial differential equations",
in: Advances in Computer Methods for Partial Differential Equations VII,
Proceedings of the 7th IMACS International Conference on Computer Methods 
for Partial Differential Equations,
Rutgers University, New Brunswick, New Jersey, June 22-24, 1992,
Eds.: R. Vichnevetsky, D. Knight and G. Richter, 
IMACS, New Brunswick, New Jersey, pp. 326-332 (1993)

                ----------------------------------

The software is also available via anonymous FTP from our ftp site:

ftp.mines.edu   or  mines.edu

FTP to that site, login with anonymous, use your email address or name
as password. Then change to the directory pub/papers/math_cs_dept/software.

The painleve software is in the subdirectory npainleve, and its 
subdirectory single.

                 ----------------------------------

In this subdirectory NPAINLEVE/SINGLE (on the disk, or on the workstation*), 
you should have the following files :

(* on the workstation there is only one directory)

In the subdirectory PROGRAM (on the disk) you will find

- np_setup.max :  the Macsyma file which initializes the program
- np_sing.max:    the Macsyma source code of the program
- np_exec.max :   the Macsyma execution file

In the subdirectory DOCUMENT (on the disk) you will find

- np_readme.txt : the file you are reading now

  and also the files for all of the papers mentioned above. 

In the subdirectory TESTDECK (on the disk) you will find

- p_kdv.com : the command file to run the Korteweg-de Vries equation
- p_kdv.dat : the data for Korteweg-de Vries equation as an example of input
- p_kdv.out : output of the Painleve test for the KdV equation

  Various other examples have been prepared. In TESTDECK
  you will find data files, command files and output files for, amongst others:

- the sine-Gordon equation (PDE): p_sg.com, p_sg.dat and p_sg.out
- the Fisher equation (in travelling frame, therefore ODE form): 
                             p_fisher.com, p_fisher.dat and p_fisher.out
- the FitzHugh-Nagumo equation (in travelling frame, thus ODE): 
                             p_fhn.dat, p_fhn.com and p_fhn.out
- the cylindrical KdV equation (PDE with variable coefficient):
  p_cylkdv.com, p_cylkdv.dat, p_cylkdv.out
- a PDE due to Sasha Mikhailov (to demonstrate the Kruskal simplification):
  p_sasha.com, p_sasha.dat and p_sasha.out

There is also a long file, called p_vareqs.dat, with various data of other
PDEs and ODEs, that were used in testing and debugging of the program.

The files np_setup.max, np_exec.max and np_sing.max, need to be copied in a 
directory on the mainframe, from where MACSYMA can batch these files.

             -----------------------------------------

How to use the program?

After the files are in place on your system (that should have Macsyma),
you may want to try the simple cases (Korteweg-de Vries, the Fisher, 
the cylindrical KdV, the FitzHugh-Nagumo, and sine-Gordon equations).
The command files and data files are all available. Also, for comparison
the output files are given. 

Once you started up Macsyma and the line (c1) comes up, just type

(c1) batch("p_kdv.com")  /* to run the KdV case */

              ---------------------------------------

                  The new features of the program.

The user can still supply a value of alpha, say -2, in the data file.
The name for alpha in the data file must be beta, to avoid confusion 
with the build in alpha. Add the following line to the data file:

beta : -2; 

-------
 
The execution of the computation of the compatibility conditions at the
resonances is still optional. If you do not want it, add the following 
line to the data file:

do_resonances: false;

If all the compatibility conditions should be computed the data file
could have the OPTIONAL line

do_resonances: true;   /* this is the default anyway */

-------

It is now possible to give uzero from outside, if so,  
the flag giving_uzero must be set to true in the data file:

giving_uzero: true; 
depends(g,x);
supplied_uzero: -12*diff(g,x)^2;

If uzero is not given the data file could have the OPTIONAL line:

giving_uzero: false;  /* this is the default anyway */

-------

Sometimes solutions of u[0] can be very messy, for those cases if you still
want to solve for u[0] explicitly set the flag solving_uzero to true
in the data file:  

solving_uzero: true;

If the program may not solve for u[0], the data file could have 
the OPTIONAL line:

solving_uzero: false;  /* this is the default anyway */

---------------

At the end of the computations several pieces of information are available:

The values of the varies coefficients in the Laurent series can be extracted 
for further manipulation as follows:
 
type uval[j,k,l], where j stands for j_th alpha, k stands for u[k], and l 
stands for l_th solution set for u[0].

The values for alpha are also available: 

type alpha[j], where j stands for j_th alpha.

The compatibility conditions are also accessible:

type compcond[j,k], where j stands for j_th alpha, k stands for k_th solution 
set for u[0].

The resonances are available:
 
type res[j,k], where j stands for j_th alpha, k stands for k_th solution 
set for u[0].

              ---------------------------------------

Note for PC-Macsyma users:

If you are using the program on a PC equipped with PC-Macsyma
you will have to specify the path where the various files are located.
For instance, if you have the files np_sing.max, np_setup.max and
np_exec.max all in the directory c:\macsyma\npainleve then 
the np_setup.max file should read:

/* Setup file for the Painleve program */

/* ********************************************************************** */
/*                           Batch file NP_SETUP.MAX                      */
/* ********************************************************************** */

  kill (all)$
  loadprint : false$
  batchload("c:\\macsyma\\npainleve\\np_sing.max")$
  derivabbrev : true$
  depends([f,eq],[t,x,y,z])$
  fx[k](x):=diff(f,x,k)$
  ftx[k,l](t,x):=diff(diff(f,t,k),x,l)$
  ftxy[k,l,m](t,x,y):=diff(diff(diff(f,t,k),x,l),y,m)$
  ftxyz[k,l,m,n](t,x,y,z):=diff(diff(diff(diff(f,t,k),x,l),y,m),z,n)$

/* *************************** END of P_SETUP.MAX *********************** */

Accordingly the batch file to run e.g. the FHN case should have the correct
path specifications for files to be batched, say from c:\macsyma\npainleve. 
For example, the command file for the FHN case should read:

/* Command file for FitzHugh-Nagumo equation in ODE form */

/* ************************************************************************ */
/*                          Batch file P_FHN.COM                            */
/* ************************************************************************ */

batch("c:\\macsyma\\npainleve\\np_setup.max")$
batch("c:\\macsyma\\npainleve\\p_fhn.dat")$
writefile("c:\\macsyma\\npainleve\\p_fhn.out");
batch("c:\\macsyma\\npainleve\\np_exec.max")$
closefile();
quit();

/* ************************** END of P_FHN.COM *************************** */

Similar changes may have to be made for other platforms. 
Consult the Macsyma guide (or manual) for your system. 

To run the FHN case at the prompt (c1) under Macsyma on your PC, type 

batch("c:\\macsyma\\npainleve\\p_fhn.com");

                 ----------------------------------

How to run examples?

After you start Macsyma on your system, type 

(c1) batch("p_kdv.com");  /* to perform the Painleve test for the KdV eq. */

This is the only line you need to type, ignoring the comment of course. 

What happens next? In turn, the batch file p_kdv.com invokes the following:

batch("np_setup.max") $        resets and initializes, before every new example
batch("p_kdv.dat") $           reads in the data file p_kdv.dat
writefile("p_kdv.out") $       opens the transcript file p_kdv.out
batch("np_exec.max") $         carries out the Painleve test for the equation
closefile() $                  closes output file p_kdv.out
quit()$                        quit from the program Macsyma

Instead of typing these lines in interactive mode you could write little
batch file (these are the xxx.com files on the disk).

           ---------------------------------------------------

Preparing the data file(s):

The terms of the equation to be tested need to be entered in a specific way:
A typical term reads ftxyz[k,l,m,n](t,x,y,z), where k, l, m, and n stand for 
the order of derivation with respect to the variables t, x, y, and z.
For ODEs the variable x is mandatory. So, one has to use : fx[.](x), where
within the brackets the order or derivation is inserted.
For simplicity, the function without derivatives may be denoted by f itself.

For example, to test the KdV one would prepare the batch file p_kdv.dat 
with one (or two lines) :

eq : ftx[1,0](t,x) + bb*f*ftx[0,1](t,x) + ftx[0,3](t,x) ;
/* where bb is an arbitrary constant */ 

The output for this ubiquitous equation is given in p_kdv.out.

Almost any letters of the alphabet can be used for arbitrary constants, 
e.g. a, b, c, d, e, ..., a1, a2, ..., b1, b2, etc. 
But there are exceptions: DO NOT USE c1, c2, ... d1, d2, ..., and e1, e2,... 
since these refer to the labels of command lines, output lines and 
intermediate expressions in Macsyma.

An exhaustive list of examples of equations are given in the file
p_vareqs.dat provided herewith. These and other examples
were used in building, testing and debugging the program.

                      ------------------------------------

A few notes about this update:

The new version 2.5 of this Painleve program allows one to use 
the simplification suggested by Martin Kruskal.

To use this option you must specify two parameters in the data file: 

do_simplification : true $
prefer_variable : x $

You may select the variables x, y, z, and t.
If indeed the variable x is preferred then the expansion variable g 
will be g(t,x,y,z) = x - h(t,y,z).

                                     dg
In the program we set automatically  -- = 1 by using gradef(g,x,1),
                                     dx
                                    
furthermore, we set gradef(g,t, - diff(h,t), etc.

The coefficients u[k] in the expansion of will be u(k,t,y,z).
For other choices of the parameter 'prefer_variable', the formula will be 
adjusted. For instance, if y were selected then g(t,x,y,z) = y - h(t,x,z).

Using this simplification drastically cuts the work. Sometimes the calculation
time is reduced to one tenth of the time it takes to carry out the original
Painleve test.

An example due to Sasha Mikhailov, that shows Kruskal's simplification,
is given in the file p_sasha.out. The interactive file that produces
this output is called p_sasha.com. The data are given in p_sasha.dat.

The fairly complicated equation comes from Jimbo-Miwa's hierarchy of evolution 
equations for the Kadomtsev Petviashvili equation and passes the Painleve test.
The traditional test (due to Weiss and others) requires about 1 hour of 
CPU time. With Kruskal's simplification the test in finished in minutes.

Note: If you are using a former version of this Painleve program already, 
then before using the new version of the program make sure you rename or 
delete the old files main_p.412 (or main_p.309), and painl412.lsp 
(or painl309.lsp) otherwise Macsyma (or you) may get confused.

          -------------------------------------

The versions of the program we provided here is for Macsyma 412 or higher
version, including Macsyma 417.125 and 2.0 for PC.

This version of the program has been tested with Macsyma version 412.61
on the VAX 11/780  of the Center for the Mathematical Sciences at the
University of Wisconsin  at Madison, Wisconsin 53705, USA and also
on the VAX 8600 of Colorado School of Mines in Golden, CO 80401, USA.
Later, tests were performed under PC-Macsyma version 417.125 on VAX 8500 and
IBM Risc 6000, and the newest version 2.0 for PC under Windows.

You need to contact me if you still have Macsyma 309 on your computer.
I can provide you with an older version of the program that is written for 
any computer that has the possibility to implement Franz Lisp MACSYMA 309.

                --------------------------------------

To learn about new updates of the program, or in case of trouble, contact me.

By phone: (303) 273-3881 (office, with messages) or 3860 (secretary), 
          or (303) 440-6089 (home, with answering service);

By fax: (303) 273-3875 (mention for Dr. Hereman)

By email: WHEREMAN@FLINT.MINES.EDU 
          or WHEREMAN@LIE.MINES.EDU

By mail:

Dr. Willy Hereman
Associate Professor
Department of Mathematical and Computer Sciences
Colorado School of Mines
Golden, Colorado 80401-1887, USA

Any comments about the program are welcomed by the authors. 
Good luck!

Willy 

Golden, August 28, 1995.

/* ************************* END of NP_README.TXT ********************** */

