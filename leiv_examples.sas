/*
leiv_examples.sas
david leonard
20 nov 2018
examples of linear errors in variables estimation
*/

/* Background
Solutions of the bivariate, linear errors-in-variables estimation problem with
unspecified errors are expected to be invariant under interchange and scaling
of the coordinates. The appealing model of normally distributed true values
and errors is unidentified without additional information.
*/

/* Purpose
The following program illustrates errors-in-variables estimation using the leiv
SAS macro. The leiv SAS macro computes the leiv posterior density of the slope.
It incorporates the fact that the slope and variance parameters together
determine the covariance matrix of the unobserved true values but is otherwise
diffuse. It is invariant to interchange and scaling of the coordinates and
depends on the data only through the sample correlation coefficient and ratio
of standard deviations.

The program illustrates the graphic and output options of the macro as well as
estimation from raw data, from summary (TYPE=CORR) data, or from sufficicient
statistics.

The program will write the results to an RTF file with the same name as the
calling program (but an RTF extension) in the current working directory.
*/

/* Reference
Leonard, David. Estimating a bivariate linear relationship.
Bayesian Anal. 6 (2011), no. 4, 727--754. doi:10.1214/11-BA627.
https://projecteuclid.org/euclid.ba/1339616542
*/


/* Current working directory and file information */
%let execfilepath=%sysget(SAS_EXECFILEPATH); /* path to the calling program */
%let execfilename=%sysget(SAS_EXECFILENAME); /* name of the calling program */
%let pwd=%substr(&execfilepath,1,%index(&execfilepath,&execfilename)-2); /* path to current working directory */
%let rtfpath=%scan(&execfilepath,1,'.').rtf; /* complete path to RTF output file */

filename curdir "&pwd"; /* specify the current working directory as an aggregate storage location */
%include curdir(leiv_macros); /* execute the macro definitions (assumes leiv_macros.sas is located in the current working directory) */

/* Examples */

options dtreset;
ods graphics / reset=all outputfmt=png noborder;
ods rtf file="&rtfpath" sasdate nogtitle nogfootnote;

title 'Examples of the LEIV SAS Macro';

title2 'Simulated Data';

%let n=20;
%let xMean=5;
%let xSD=4;
%let exSD=5;
%let eySD=3;
%let beta=1;
%let beta0=2;
%let seed=1123;

title3 "N=&n, Intercept=&beta0, Slope=&beta";
data test (keep=id xObserved yObserved);
do id=1 to &n;
	xTrue=&xMean+&xSD*rannor(&seed);
	yTrue=&beta0+&beta*xTrue;
	xObserved=xTrue+&exSD*rannor(0);
	yObserved=yTrue+&eySD*rannor(0);
	output;
end;
label
xObserved='Observed value of x'
yObserved='Observed value of y'
;
run;
proc print data=test noobs;
run;

title4 'Estimation directly from the data';
%leiv(x=xObserved,y=yObserved,data=test,suffdata=suffstats,estdata=estimates,plots=scatter density,postdata=posterior);

title4 'Estimation from the associated TYPE=CORR summary data';
proc corr data=test noprint out=testcorr nomiss;
var xObserved yObserved;
run;
%leiv(x=xObserved,y=yObserved,data=testcorr);

title2 'Calculate a LEIV posterior density from sufficient statistics';
%let n=20;
%let sdRatio=0.68;
%let cor=0.70;
data dummycorr (type=corr);
_TYPE_='MEAN'; _NAME_=''; x=0; y=0;
output;
_TYPE_='STD';_NAME_=''; x=1; y=&sdRatio;
output;
_TYPE_='N'; _NAME_=''; x=&n; y=&n;
output;
_TYPE_='CORR'; _NAME_='x'; x=1; y=&cor;
output;
_TYPE_='CORR'; _NAME_='y'; x=&cor; y=1;
output;
run;
%leiv(x=x,y=y,data=dummycorr,suffdata=dummysuff,estdata=dummyest,plots=density,postdata=dummypost,postlower=-1,postupper=3,poststep=0.02,level=0.9);

title;

ods rtf close;
