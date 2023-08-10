%-----------------------------------------------------------------------------------------------------------------------------------
% solves for the slip profile and outlet pressure needed for a desired deflection profile for "Thin substrate with given platform 
% profile and given axial variation of Lame parameters"
%-----------------------------------------------------------------------------------------------------------------------------------
function tapersoln = tapersolve(inversesystem,simul)

%-----------------------------------------------------------------------------------------------------------------------------------
%	loading system variables
	struct2vars(inversesystem);
%	loading simulation variables
	struct2vars(simul);
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
	% system variables
	L =				L1;
	bet =			Delt/L;
	gamm =			H/L;
	kapp =			bslipper/L;
	ph0 =			((bet*kapp)/(gamm^3))*((muvisc*Q)/((min(lambd+2*G))*L^2));
	tildph =		zeros(1,nx);
	for ix = 1:nx
		tildph(ix)=	min(lambd+2*G)/(lambd(ix)+2*G(ix));
	end
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	obtaining slip profile
	xigreek =		zeros(1,nx);
	Abarph =		zeros(nx,nx);
	bbarph =		zeros(1,nx);
	ix = 1;
	atemp1 =							a('x','fd',ix,ix,x,x);
	atemp2 =							a('','fd',ix,ix,x,x);
	for jx = 0:2
		Abarph(ix,ix+jx) =				atemp1(jx+3,0+3)*h(ix)+atemp2(jx+3,0+3)*differential1d(h,'x','fd',ix,ix,'x',x,x);
	end
	bbarph(ix) =						(12.0*(1.0-((ph0*h(ix))/gamm)+((bslip(ix))/(gamm*L))))/ ...
										((1.0-((ph0*h(ix))/gamm)+((4*bslip(ix))/(gamm*L)))*((1.0-((ph0*h(ix))/gamm))^3));
	for ix = 2:nx-1
		atemp1 =						a('x','cd',ix,ix,x,x);
		atemp2 =						a('','cd',ix,ix,x,x);
		for jx = -1:1
			Abarph(ix,ix+jx) =			atemp1(jx+3,0+3)*h(ix)+atemp2(jx+3,0+3)*differential1d(h,'x','cd',ix,ix,'x',x,x);
		end
		bbarph(ix) =					(12.0*(1.0-((ph0*h(ix))/gamm)+((bslip(ix))/(gamm*L))))/ ...
										((1.0-((ph0*h(ix))/gamm)+((4*bslip(ix))/(gamm*L)))*((1.0-((ph0*h(ix))/gamm))^3));
	end
	ix = nx;
	Abarph(ix,ix) =						1.0;
	bbarph(ix) =						1.0;
	barph =								(Abarph\bbarph')';
	for ix = 1:nx
		xigreek(ix) =					tildph(ix)/barph(ix);
	end
	if (min(xigreek) < 0)
		xigreek =						xigreek-min(xigreek);
	end
	scalexi =							1/max(xigreek);
	xigreek =							xigreek*scalexi;
	Delt =								Delt/scalexi;
	p0 =								-(kapp/(gamm^3))*((muvisc*Q)/(L^2))*(barph(nx)*h(nx));
%-----------------------------------------------------------------------------------------------------------------------------------
	
%-----------------------------------------------------------------------------------------------------------------------------------
%	feeding variables to output struct
	tapersoln.xigreek =	xigreek;
	tapersoln.Delt =	Delt;
	tapersoln.p0 =		p0;
%-----------------------------------------------------------------------------------------------------------------------------------
end