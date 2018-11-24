/*
leiv_macros.sas
david leonard
20 nov 2018
linear errors in variables estimation
*/

/* Background
Solutions of the bivariate, linear errors-in-variables estimation problem with
unspecified errors are expected to be invariant under interchange and scaling
of the coordinates. The appealing model of normally distributed true values
and errors is unidentified without additional information.
*/

/* Purpose
The following program provides SAS macros for computing the leiv posterior
density of the slope. It incorporates the fact that the slope and variance
parameters together determine the covariance matrix of the unobserved true
values but is otherwise diffuse. It is invariant to interchange and scaling of
the coordinates and depends on the data only through the sample correlation
coefficient and ratio of standard deviations.

Graphic and output options are available as are methods for estimation from
raw data, from summary (TYPE=CORR) data, or from sufficient statistics.
*/

/* Reference
Leonard, David. Estimating a bivariate linear relationship.
Bayesian Anal. 6 (2011), no. 4, 727--754. doi:10.1214/11-BA627.
https://projecteuclid.org/euclid.ba/1339616542
*/


/* Utility macros for Note, Warning and Error logging */
%macro message(msg);
	%put &msg.%str(.);
	%put;
%mend message;
%macro note(msg);
	%message(NOTE: &msg.);
%mend note;
%macro warning(msg);
	%message(WARNING: &msg.);
%mend warning;
%macro error(msg);
	%message(ERROR: &msg.);
	%abort;
%mend error;


/* Constants used in LEIV estimation */
%let pi_2=%sysevalf(2*%sysfunc(atan(1))); /* pi/2 */


/* Main leiv macro */
%macro leiv(y=,x=,data=,suffdata=,estdata=,plots=,postdata=,postlower=,postupper=,poststep=,level=,refresh=,reltol=,abstol=,jmax=,k=);

/*
arguments
---------------------------------------------------------------------------------------------------------------------------
y			|	(required) identifies the y variable
x			|	(required) identifies the x variable
data		|	(required) identifies a SAS dataset including the x and y variables or the associated SAS TYPE=CORR dataset 
suffdata	|	(optional) names an output dataset containing the sufficient statistics
estdata		|	(optional) names an output dataset containing the posterior estimates
plots		|	(optional) a space separated list with possible entries: none scatter density (none by default)
postdata	|	(optional) names an output dataset containing a uniform sample of the posterior density
postlower	|	(optional) the lower limit for sampling the posterior density
postupper	|	(optional) the upper limit for sampling the posterior density
poststep	|	(optional) the step size for sampling the posterior density
level		|	(optional) the probability level of the shortest probability level (0.95 by default)
refresh		|	(optional) if 1 then refreshes LEIV functions defined previously (0 by default)
reltol      |   (optional) relative tolerance for numerical integration using Romberg's method [Numerical Recipes 4.3]
abstol      |   (optional) absolute tolerance for numerical integration using Romberg's method [Numerical Recipes 4.3]
jmax        |   (optional) 2^(jmax-1) is the maximum allowed number of steps in Romberg's method [Numerical Recipes 4.3]
k           |   (optional) number of points in the polynomial extrapolation in Romberg's method [Numerical Recipes 4.3]
---------------------------------------------------------------------------------------------------------------------------
*/

/* interpret the call */
%let suffdataout=%length(&suffdata) gt 0;
%if not &suffdataout %then
	%let suffdata=sufftmp_;
%let estdataout=%length(&estdata) gt 0;
%if not &estdataout %then
	%let estdata=esttmp_;
%let scatterplot_=0;
%let densityplot_=0;
%let item=1;
%let plotitem=%scan(&plots,&item,%str( ));
%do %while(%length(&plotitem) gt 0);
	%if &plotitem eq scatter %then
		%let scatterplot_=1;
	%else %if &plotitem eq density %then
		%let densityplot_=1;
	%else %if &plotitem ne none %then
		%warning(Unknown plot request &=plotitem%str(;) expecting NONE%str(,) SCATTER or DENSITY);
	%let item=%eval(&item+1);
	%let plotitem=%scan(&plots,&item,%str( ));
%end;
%let postdataout=%length(&postdata) gt 0;
%let posterior=%eval(&postdataout or &densityplot_);
%if &posterior and not &postdataout %then
	%let postdata=posttmp_;
%if %length(&postlower) eq 0 %then
	%let postlower=.;
%if %length(&postupper) eq 0 %then
	%let postupper=.;
%if %length(&poststep) eq 0 %then
	%let poststep=.;
%if %length(&level) eq 0 %then
	%let level=0.95;
%if %length(&refresh) eq 0 %then
	%let refresh=0;
%if %length(&reltol) eq 0 %then
	%let reltol=%sysevalf(%sysfunc(CONSTANT(SQRTMACEPS))**0.5);
%if %length(&abstol) eq 0 %then
	%let abstol=&reltol;
%if %length(&jmax) eq 0 %then
	%let jmax=20;
%if %length(&k) eq 0 %then
	%let k=5;

/* define LEIV functions */
%let leivlib=%sysfunc(exist(leivfuncs));
%if &leivlib and not &refresh %then
	%note(Using pre-existing LEIV functions%str(;) use the LEIV option REFRESH=1 to refresh);
%else %do;
	%if &leivlib %then %do;
		%note(Deleting LEIV functions defined previously);
		proc delete data=leivfuncs;
		run;
	%end;
	%note(Defining LEIV functions);
	options cmplib=work.leivfuncs;
	proc fcmp outlib=work.leivfuncs.LEIV;

	subroutine polint(xary[*],yary[*],n,x,y,dy,rc); /* utility routine called by Romberg integrators but not meant to be called by the user */
	/* return the interpolating value, y, at x and an error estimate, dy of the polynomial of degree n-1 through the n given points */
		outargs y, dy, rc, xary, yary; /* pass the x and y arrays by reference for efficiency */
		
		array c[1] / nosym;
		call dynamic_array(c,n);
		array d[1] / nosym;
		call dynamic_array(d,n);

		rc = 0; /* assume success */
		dif=abs(x-xary[1]);
		ns=1;
		do i=1 to n;
			dift=abs(x-xary[i]);
			if dift lt dif then do;
				ns=i;
				dif=dift;
			end;
			c[i]=yary[i];
			d[i]=yary[i];
		end;

		y=yary[ns]; /* initial approximation to y */
		ns + -1;
		do m=1 to n-1;
			do i=1 to n-m;
				ho=xary[i]-x;
				hp=xary[i+m]-x;
				w=c[i+1]-d[i];
				den=ho-hp;
				if den ne 0 then do;
					den=w/den;
					d[i]=hp*den;
					c[i]=ho*den;
				end;
				else
					rc=1; /* failure: 2 x's are equal to within roundoff */
			end;
			if 2*ns lt n-m then
				dy=c[ns+1];
			else do;
				dy=d[ns];
				ns + -1;
			end;
			y+dy;
		end;
	endsub;

	function tIntegrand(t,tLower,tUpper,v);
	*	tBig=squantile('T',&abstol,v)**(1-0.5**v); /* accurate but slow */
		tBig=exp(30/v**0.8)**(1-0.5**v);
		if tUpper le tBig then do;
			if t lt tUpper then do;
				vm=v-1;
				vp=v+1;
				F=vm/vp*(v+t*t)/(tUpper+t-2*tLower)/(tUpper-t);
				if F gt 0 then
					ft=pdf('T',t,v)*cdf('F',F,vp,vm);
				else
					ft=0;
				return(ft);
			end;
			else if t eq tUpper then do;
				ft=pdf('T',t,v);
				return(ft);
			end;
		end;
		else
			return(0);
	endsub;

	function trapt(lower,upper,n,tLower,tUpper,v); /* utility function called by qrombt but not meant to be called by the user */
	/* returns the nth stage of refinement of an extended trapezoid rule integration of the specified function */
		static s;
		if n eq 1 then
			s=0.5*(upper-lower)*(tIntegrand(lower,tLower,tUpper,v)+tIntegrand(upper,tLower,tUpper,v));
		else do;
			j=2**(n-2);
			del=(upper-lower)/j;
			x=lower+0.5*del;
			sum=0;
			do i=1 to j;
				sum+tIntegrand(x,tLower,tUpper,v);
				x+del;
			end;
			s=0.5*(s+(upper-lower)*sum/j);
		end;
		return(s);
	endsub;

	function qrombt(lower,upper,tLower,tUpper,v,rc); /* utility function used to calculate the central integral but not meant to be called by the user */
	/* integrate a function over a closed interval using Romberg's method
	of order 2K, where, for example, K=2 is simpson's rule */
	outargs rc;

		array s[1] / nosym;
		call dynamic_array(s,&jmax); /* successive trapezoidal approximations */
		array h[1] / nosym;
		call dynamic_array(h,&jmax+1); /* successive relative stepsizes */

		array sk[&k] / nosym;
		array hk[&k] / nosym;

		h[1]=1.0;
		j=1;
		rc=1; /* assume failure */
		do while(j le &jmax and rc ne 0);
			s[j]=trapt(lower,upper,j,tLower,tUpper,v);
			if j ge &k then do;
				do i=1 to &k;
					hk[i]=h[j-&k+i];
					sk[i]=s[j-&k+i];
				end;
				result=.;
				dres=.;
				rcode=.;
				call polint(hk,sk,&k,0.0,result,dres,rcode);
				if rcode eq 0 then do;
					abserr=abs(dres);
					if abserr le &reltol*abs(result) and abserr le &abstol then
						rc=0; /* success */
				end;
			end;
			h[j+1]=0.25*h[j];
			j+1;
		end;

		return(result);
	endsub;

	function I(b,r,s,v);
	/* central integral of dimensionless, scalar arguments */
		tLower=-r/s;
		tUpper=(b-r)/s;
		rcode=.;
		result=qrombt(tLower,tUpper,tLower,tUpper,v,rcode);
		if rcode eq 0 then
			return(result);
		else
			return(.);
	endsub;

	function J(b,cor,s,v);
	/* reduced likelihood */
		if b ne 0 then do;
			bSign=sign(b);
			bAbs=bSign*b;
			rbSign=cor*bSign;
			return(I(bAbs,rbSign,s,v)+I(1/bAbs,rbSign,s,v));
		end;
		else
			return(0);
	endsub;

	function p0(b);
	/* cauchy prior density */
		return(pdf('CAUCHY',b));
	endsub;

	function p1(b,cor,s,v);
	/* unnormalized posterior density */
		return(J(b,cor,s,v)*p0(b));
	endsub;

	function qIntegrand(q,cor,s,v);
	/* unnormalized posterior density in terms of q=atan(b) */
		if -&pi_2 lt q lt &pi_2 then do;
			b=tan(q);
	*		Jac=sec(q)**2; /* Jacobian |db/dq| */
	*		return(p1(b,cor,s,v)*Jac);
			return(J(b,cor,s,v)/&pi_2/2); /* p0(b)*Jac = 1/pi */
		end;
		else
			return(0);
	endsub;

	function trapq(lower,upper,n,cor,s,v); /* utility function called by qrombq but not meant to be called by the user */
	/* returns the nth stage of refinement of an extended trapezoid rule integration of the specified function */
		static r;
		if n eq 1 then
			r=0.5*(upper-lower)*(qIntegrand(lower,cor,s,v)+qIntegrand(upper,cor,s,v));
		else do;
			j=2**(n-2);
			del=(upper-lower)/j;
			x=lower+0.5*del;
			sum=0;
			do i=1 to j;
				sum+qIntegrand(x,cor,s,v);
				x+del;
			end;
			r=0.5*(r+(upper-lower)*sum/j);
		end;
		return(r);
	endsub;

	function qrombq(lower,upper,cor,s,v,rc); /* utility function used to calculate the cumulative distribution but not meant to be called by the user */
	/* integrate a function over a closed interval using Romberg's method
	of order 2K, where, for example, K=2 is simpson's rule */
	outargs rc;

		array r[1] / nosym;
		call dynamic_array(r,&jmax); /* successive trapezoidal approximations */
		array h[1] / nosym;
		call dynamic_array(h,&jmax+1); /* successive relative stepsizes */

		array rk[&k] / nosym;
		array hk[&k] / nosym;

		h[1]=1.0;
		j=1;
		rc=1; /* assume failure */
		do while(j le &jmax and rc ne 0);
			r[j]=trapq(lower,upper,j,cor,s,v);
			if j ge &k then do;
				do i=1 to &k;
					hk[i]=h[j-&k+i];
					rk[i]=r[j-&k+i];
				end;
				result=.;
				dres=.;
				rcode=.;
				call polint(hk,rk,&k,0.0,result,dres,rcode);
				if rcode eq 0 then do;
					abserr=abs(dres);
					if abserr le &reltol*abs(result) and abserr le &abstol then
						rc=0; /* success */
				end;
			end;
			h[j+1]=0.25*h[j];
			j+1;
		end;

		return(result);
	endsub;

	function cdfb(b,cor,s,v,nconst);
	/* normalized posterior cumulative distribution function */
		q=atan(b);
		lower=-&pi_2;
		rcode=.;
		cdf=qrombq(lower,q,cor,s,v,rcode)/nconst;
		if rcode eq 0 then
			return(cdf);
		else
			return(.);
	endsub;

	function p50(cor,s,v,nconst,rc);
	/* posterior median */
		outargs rc;
		array opts[5] / nosym;
		opts[1]=cor; /* initial value */
		opts[2]=&abstol; /* absolute criterion */
		opts[3]=&reltol; /* relative criterion */
		opts[4]=100; /* maximum number of iterations */
		opts[5]=.; /* return code: 0=success, 1=doesn't converge, 2=can't calculate, 3=maximum iterations exceeded, 4=initial objective missing */
		rcode=.;
		result=solve('cdfb',opts,0.5,.,cor,s,v,nconst);
		rc=opts[5];
		return(result);
	endsub;

	function upper(lower,median,cor,s,v); /* utility function called by prob and probInt but not meant to be called by the user */
	/* upper limit of the shortest probability interval as a function of the lower limit */
		if cor gt 0 then
			init=1/cor+2/sqrt(v)*(1/cor-median);
		else if cor lt 0 then
			init=cor+2/sqrt(v)*(cor-median);
		else
			init=1+2/sqrt(v)*(1-median);
		array opts[5] / nosym;
		opts[1]=init; /* initial value */
		opts[2]=&abstol; /* absolute criterion */
		opts[3]=&reltol; /* relative criterion */
		opts[4]=100; /* maximum number of iterations */
		opts[5]=.; /* return code: 0=success, 1=doesn't converge, 2=can't calculate, 3=maximum iterations exceeded, 4=initial objective missing */
		b=solve('p1',opts,p1(lower,cor,s,v),.,cor,s,v);
		return(b);
	endsub;

	function prob(lower,median,cor,s,v,nconst); /* utility function called by probInt but not meant to be called by the user */
	/* probability between the upper and lower limits of a shortest probability interval */
		upper=upper(lower,median,cor,s,v);
		return(cdfb(upper,cor,s,v,nconst)-cdfb(lower,cor,s,v,nconst));
	endsub;

	subroutine probInt(level,lower,upper,median,cor,s,v,nconst,rc);
	/* shortest (100*level)% probability interval */
		outargs lower, upper, rc;
		if cor gt 0 then
			b=cor;
		else if cor lt 0 then
			b=1/cor;
		else
			b=1;
		width=2/sqrt(v)*(median-b);
		lower=b-width;
		p12=prob(lower,median,cor,s,v,nconst);
		do while(2*p12 lt 1+level); /* initialize with an interval wider than desired */
			width=2*width;
			lower=b-width;
			p12=prob(lower,median,cor,s,v,nconst);
		end;
		array opts[5] / nosym;
		opts[1]=lower; /* initial value */
		opts[2]=&abstol; /* absolute criterion */
		opts[3]=&reltol; /* relative criterion */
		opts[4]=100; /* maximum number of iterations */
		opts[5]=.; /* return code: 0=success, 1=doesn't converge, 2=can't calculate, 3=maximum iterations exceeded, 4=initial objective missing */
		lower=solve('prob',opts,level,.,median,cor,s,v,nconst);
		upper=upper(lower,median,cor,s,v);
		rc=opts[5];
	endsub;

	run;

	%note(Ending LEIV function definitions);
%end;

/* LEIV titles */
%let leivtitlefont="sans-serif" bold italic;
%let leivwarningfont="sans-serif";
proc sql noprint;/* get the current title number */
select max(number)
	into :TitleCurrent
from Dictionary.Titles
where type eq 'T'; /* exclude footnotes (type='F') */
quit;
title%eval(&TitleCurrent+2) font=&leivtitlefont 'Linear Errors-In-Variables Estimation';

%note(Beginning LEIV estimation);
/* calculate and print the sufficient statistics */
%if %length(&x) eq 0 %then
	%error(Missing &=x);
%if %length(&y) eq 0 %then
	%error(Missing &=y);
%let dsid = %sysfunc(open(&data));
%let type_=%sysfunc(attrc(&dsid,TYPE));
%let rc = %sysfunc(close(&dsid));
%if %length(&type_) eq 0 %then %do;
	/* check for valid data */
	%let dsid = %sysfunc(open(&data));
	%let n_=%sysfunc(attrn(&dsid,NOBS));
	%let rc = %sysfunc(close(&dsid));
	%if &n_ lt 2 %then
		%error(Requires n>=2 data points);
	data datatmp_;
	set &data;
	if not missing(&x) and not missing(&y);
	run;
	%let dsid = %sysfunc(open(datatmp_));
	%let n_=%sysfunc(attrn(&dsid,NOBS));
	%let rc = %sysfunc(close(&dsid));
	%if &n_ lt 2 %then
		%error(Requires n>=2 non-missing data points);
	/* calculate summary statistics */
	proc sql noprint;
	create table &suffdata as
	select
		sum(1) as n,
		sum(&x)/(calculated n) as xMean,
		sum(&y)/(calculated n) as yMean,
		sum(&x*&x)/(calculated n)-(calculated xMean)**2 as Sxx,
		sum(&y*&y)/(calculated n)-(calculated yMean)**2 as Syy,
		sum(&x*&y)/(calculated n)-(calculated xMean)*(calculated yMean) as Sxy
	from datatmp_;
	quit;
	data &suffdata;
	set &suffdata;
	call symput('Sxx_',Sxx);
	call symput('Syy_',Syy);
	run;
	%if %sysevalf(&Sxx_ eq 0 and &Syy_ eq 0, boolean) %then
		%error(Requires n>=2 distinct data points);
	%if %sysevalf(&Syy_ eq 0, boolean) %then
		%error(Singular data%str(;) linear with zero slope);
	%if %sysevalf(&Sxx_ eq 0, boolean) %then
		%error(Singular data%str(;) linear with infinite slope);
	/* calculate sufficient statistics */
	data &suffdata (keep=n xMean yMean sdRatio cor);
	set &suffdata;
	sdRatio=sqrt(Syy/Sxx);
	cor=Sxy/sqrt(Sxx*Syy);
	call symput('cor_',cor);
	run;
%end;
%else %if &type_ eq CORR %then %do;
	/* read summary statistics */
	data &suffdata;
	set &data (where=(_type_ eq 'N'));
	nx=&x;
	ny=&y;
	call symput('nx_',nx);
	call symput('ny_',ny);
	run;
	%if %sysevalf(&nx_ ne &ny_, boolean) %then
		%error(Unequal sample size in &x and &y%str(;) use the NOMISS option in the PROC CORR statement);
	%if %sysevalf(&nx_ lt 2, boolean) %then
		%error(Requires non-missing n>=2);
	data &suffdata;
	set &suffdata;
	set &data (where=(_type_ eq 'MEAN'));
	xMean=&x;
	yMean=&y;
	set &data (where=(_type_ eq 'STD'));
	SDx=&x;
	SDy=&y;
	call symput('SDx_',SDx);
	call symput('SDy_',SDy);
	run;
	%if %sysevalf(&SDx_ lt 0 or &SDy_ lt 0, boolean) %then
		%error(Requires standard deviations >= 0);
	%if %sysevalf(&SDx_ eq 0 and &SDy_ eq 0, boolean) %then
		%error(Requires n>=2 distinct data points);
	%if %sysevalf(&SDy_ eq 0, boolean) %then
		%error(Singular data%str(;) linear with zero slope);
	%if %sysevalf(&SDx_ eq 0, boolean) %then
		%error(Singular data%str(;) linear with infinite slope);
	/* calculate sufficient statistics */
	data &suffdata (keep=n xMean yMean sdRatio cor);
	set &suffdata;
	n=floor(nx);
	call symput('n_',n);
	sdRatio=SDy/SDx;
	set &data (where=(_type_ eq 'CORR' and _name_ eq "&x"));
	cor=&y;
	call symput('cor_',cor);
	run;
	%if %sysevalf(&cor_ lt -1 or &cor_ gt 1, boolean) %then
		%error(Requires correlation in [-1%str(,)1]);
	%if %sysevalf(&n_ eq 2 and (&cor_ ne -1 and &cor_ ne 1), boolean) %then
		%error(Requires cor=-1 or cor=1 if n=2);
	%if &scatterplot_ eq 1 %then %do;
		%warning(Scatterplot request is ignored for TYPE=CORR data);
		%let scatterplot_=0;
	%end;
%end;
%else
	%error(Unrecognized data type=&type_);
title%eval(&TitleCurrent+3) font=&leivtitlefont 'Sufficient Statistics';
proc print data=&suffdata noobs;
var n xMean yMean sdRatio cor;
format n 8.;
run;

%if %sysevalf(-1 lt &cor_ and &cor_ lt 1, boolean) %then %do;

	/* calculate intermediate statistics used later */
	data sufftmp1_;
	set &suffdata;
	v=n-1;
	s=sqrt((1-cor*cor)/v);
	run;

	/* calculate the normalization constant */
	data sufftmp1_ (keep=n xMean yMean sdRatio cor s v normconst);
	set sufftmp1_;
	q1=-&pi_2;
	if cor gt 0 then do;
		q2=0;
		q3=atan(cor);
		q4=atan(1/cor);
	end;
	else if cor lt 0 then do;
		q2=atan(1/cor);
		q3=atan(cor);
		q4=0;
	end;
	else do;
		q2=atan(-1);
		q3=0;
		q4=atan(1);
	end;
	q5=&pi_2;
	/* integrate over the partition (q1,q2)(q2,q3)(q3,q4)(q4,q5) */
	rcode=.;
	k12=qrombq(q1,q2,cor,s,v,rcode);
	if rcode eq 0 then do;
		k23=qrombq(q2,q3,cor,s,v,rcode);
		if rcode eq 0 then do;
			k34=qrombq(q3,q4,cor,s,v,rcode);
			if rcode eq 0 then
				k45=qrombq(q4,q5,cor,s,v,rcode);
		end;
	end;
	if rcode eq 0 then
		normconst=k12+k23+k34+k45;
	call symput('rcode_',rcode);
	run;
	%if %sysevalf(&rcode_ ne 0, boolean) %then
		%error(Posterior density normalization failed%str(;) try reducing &=abstol);

	/* calculate and print the parameter estimates */
	title%eval(&TitleCurrent+3) font=&leivtitlefont 'Posterior Estimates';
	data &estdata;
	set sufftmp1_;
	rcode_median=.;
	bmedian=p50(cor,s,v,normconst,rcode_median); /* posterior median */
	if rcode_median le 0 then do;
		Slope=bmedian*sdRatio;
		Intercept=yMean-Slope*xMean;
		call symput('Slope_',Slope);
		call symput('Intercept_',Intercept);
	end;
	call symput('rcode_median',rcode_median);
	run;
	%if %sysevalf(&rcode_median ne 0, boolean) %then
		%error(Posterior median estimation failed%str(;) return code = %kcmpres(&rcode_median));
	data &estdata (keep=intercept slope lowerPL upperPL);
	set &estdata;
	blower=.;
	bupper=.;
	rcode_interval=.;
	call probInt(&level,blower,bupper,bmedian,cor,s,v,normconst,rcode_interval);
	if rcode_interval eq 0 then do;
		LowerPL=blower*sdRatio; label LowerPL="Lower %sysevalf(100*&level)% PL";
		UpperPL=bupper*sdRatio; label UpperPL="Upper %sysevalf(100*&level)% PL";
		call symput('LowerPL_',LowerPL);
		call symput('UpperPL_',UpperPL);
	call symput('rcode_interval',rcode_interval);
	end;
*	OLSyx=cor*sdRatio; /* ordinary least squares regression of y on x */
*	OLSxy=1/(cor/sdRatio); /* reciprocal of ordinary least squares regression of x on y */
	run;
	%if %sysevalf(&rcode_median ne 0, boolean) %then
		%warning(Posterior probability interval estimation failed%str(;) return code = %kcmpres(&rcode_interval));
	proc print data=&estdata noobs label;
	var Intercept Slope LowerPL UpperPL;
	run;
	%if &scatterplot_ %then %do;
		title%eval(&TitleCurrent+3) font=&leivtitlefont 'Scatter plot with posterior median fitted line';
		proc sgplot data=datatmp_ noautolegend;
		scatter x=&x y=&y / markerattrs=(color=gray);
		lineparm x=0 y=&Intercept_ slope=&Slope_ / clip lineattrs=(color=black thickness=2) name='fit' legendlabel='Posterior Median Fit';
		*keylegend 'fit';
		xaxis labelattrs=(weight=bold);
		yaxis labelattrs=(weight=bold);
		run;
	%end;

	%if &posterior %then %do;
		/* sample the posterior density */
		data &postdata (keep=Slope Density);
		set sufftmp1_;
		if cor gt 0 then do;
			min=sdRatio*cor-20/sqrt(v)*(&Slope_-sdRatio*cor);
			max=sdRatio/cor+20/sqrt(v)*(sdRatio/cor-&Slope_);
		end;
		else if cor lt 0 then do;
			min=sdRatio/cor-20/sqrt(v)*(&Slope_-sdRatio/cor);
			max=sdRatio*cor+20/sqrt(v)*(sdRatio*cor-&Slope_);
		end;
		else do;
			min=-1-20/sqrt(v)*(&Slope_+1);
			max=1+20/sqrt(v)*(1-&Slope_);
		end;
		lower=&postlower;
		upper=&postupper;
		step=&poststep;
		if missing(lower) then
			lower=min;
		if missing(upper) then
			upper=max;
		if missing(step) then
			step=(upper-lower)/100;
		do Slope=lower to upper by step;
			Density=p1(Slope/sdRatio,cor,s,v)/normconst/sdRatio;
			output;
		end;
		label Density='Posterior Density';
		run;
		%if &densityplot_ %then %do;
			title%eval(&TitleCurrent+3) font=&leivtitlefont "Posterior density with median and shortest %sysevalf(100*&level)% probability interval";
			proc sgplot data=&postdata;
			series x=Slope y=Density / lineattrs=(color=black thickness=2);
			refline &Slope_ / axis=x lineattrs=(color=black pattern=shortdash thickness=1);
			refline &LowerPL_ / axis=x lineattrs=(color=black pattern=dot thickness=1);
			refline &UpperPL_ / axis=x lineattrs=(color=black pattern=dot thickness=1);
			xaxis labelattrs=(weight=bold);
			yaxis labelattrs=(weight=bold);
			run;
		%end;
	%end;

	proc delete data=sufftmp1_;
	run;

%end;
%else %do;
	%warning(Exactly linear data);
	title%eval(&TitleCurrent+3) font=&leivwarningfont 'Warning: Exactly linear data';

	/* calculate and print the parameter estimates */
	title%eval(&TitleCurrent+4) font=&leivtitlefont 'Posterior Estimates';
	data &estdata (keep=intercept slope lowerPL upperPL);
	set &suffdata;
	bmedian=cor;
	Slope=bmedian*sdRatio;
	Intercept=yMean-Slope*xMean;
	call symput('Slope_',Slope);
	call symput('Intercept_',Intercept);
	blower=cor;
	bupper=cor;
	LowerPL=blower*sdRatio; label LowerPL="Lower %sysevalf(100*&level)% PL";
	UpperPL=bupper*sdRatio; label UpperPL="Upper %sysevalf(100*&level)% PL";
		call symput('LowerPL_',LowerPL);
		call symput('UpperPL_',UpperPL);
	OLSyx=cor*sdRatio;
	OLSxy=1/(cor/sdRatio);
	run;
	proc print data=&estdata noobs label;
	var Intercept Slope LowerPL UpperPL;
	run;
	%if &scatterplot_ %then %do;
		title%eval(&TitleCurrent+4) font=&leivtitlefont 'Scatter plot with posterior median fitted line';
		proc sgplot data=datatmp_ noautolegend;
		scatter x=&x y=&y / markerattrs=(color=gray);
		lineparm x=0 y=&Intercept_ slope=&Slope_ / clip lineattrs=(color=black thickness=2) name='fit' legendlabel='Posterior Median Fit';
		*keylegend 'fit';
		xaxis labelattrs=(weight=bold);
		yaxis labelattrs=(weight=bold);
		run;
	%end;

	%if &posterior %then %do;
		/* sample the posterior density */
		data &postdata (keep=Slope Density);
		set &suffdata;
		SlopeEst=cor*sdRatio;
		min=SlopeEst-1;
		max=SlopeEst+1;
		lower=&postlower;
		upper=&postupper;
		step=&poststep;
		if missing(lower) then
			lower=min;
		if missing(upper) then
			upper=max;
		if missing(step) then
			step=(upper-lower)/100;
		do Slope=lower to upper by step;
			Density=Slope eq SlopeEst;
			output;
		end;
		Slope=SlopeEst;
		Density=1;
		output;
		label Density='Posterior Density';
		run;
		proc sort data=&postdata nodup;
		by Slope;
		run;
		%if &densityplot_ %then %do;
			title%eval(&TitleCurrent+4) font=&leivtitlefont 'Posterior density';
			proc sgplot data=&postdata (where=(Slope eq &Slope_));
			scatter x=Slope y=Density / markerattrs=(color=black symbol=circlefilled);
			xaxis labelattrs=(weight=bold);
			yaxis labelattrs=(weight=bold);
			run;
		%end;
	%end;


%end;

/* clean up */
%if %length(&type_) eq 0 %then %do;
	proc delete data=datatmp_;
	run;
%end;
%if not &suffdataout %then %do;
	proc delete data=sufftmp_;
	run;
%end;
%if not &estdataout %then %do;
	proc delete data=esttmp_;
	run;
%end;
%if &posterior and  not &postdataout %then %do;
	proc delete data=posttmp_;
	run;
%end;

title%eval(&TitleCurrent+1);

%note(Ending LEIV estimation);

%mend leiv;
