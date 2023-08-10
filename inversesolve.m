%-----------------------------------------------------------------------------------------------------------------------------------
% solves for the slip profile and outlet pressure needed for a desired deflection profile for "Thin substrate with given platform 
% profile and given axial variation of Lame parameters"
%-----------------------------------------------------------------------------------------------------------------------------------
function inversesoln = inversesolve(inversesystem,simul)

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
	ph0 =			((bet*kapp)/(gamm^3))*((muvisc*Q)/(min(lambd+2*G))*L^2);
	tildph =		zeros(1,nx);
	barph =			zeros(1,nx);
	for ix = 1:nx
		tildph(ix)=	min(lambd+2*G)/(lambd(ix)+2*G(ix));
		barph(ix) =	tildph(ix)/xigreek(ix);
	end
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	obtaining derivatives of desired profile
	hbarph =				zeros(1,nx);
	for ix = 1:nx
		hbarph(ix) =		h(ix)*barph(ix);
	end
	dhbarphdx =				zeros(1,nx);
	ix = 1;
	dhbarphdx(1,ix) =		differential1d(hbarph,'x','fd',ix,ix,'x',x,x);
	for ix = 2:nx-1
		dhbarphdx(1,ix) =	differential1d(hbarph,'x','cd',ix,ix,'x',x,x);
	end
	ix = nx;
	dhbarphdx(1,ix) =		differential1d(hbarph,'x','bd',ix,ix,'x',x,x);
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%	obtaining slip profile
	bslip =			zeros(1,nx);
	for ix = 1:nx
		bslip(ix) =	-(((gamm-ph0*h(ix))*L)/4)* ...
					((12*gamm^3-(gamm-ph0*h(ix))^3*dhbarphdx(ix))/(3*gamm^3-(gamm-ph0*h(ix))^3*dhbarphdx(ix)));
	end
	if (singularfix == 1)
		for i = 1:5
			for ix = 1:nx
				if (ix == 1)
					if (abs(3*gamm^3-(gamm-ph0*h(ix))^3*dhbarphdx(ix)) < 1e-1)
						bslip(ix) =		2*bslip(ix+1)-bslip(ix+2);
					end
				elseif (ix == nx)
					if (abs(3*gamm^3-(gamm-ph0*h(ix))^3*dhbarphdx(ix)) < 1e-1)
						bslip(ix) =		2*bslip(ix-1)-bslip(ix-2);
					end
				else
					if (abs(3*gamm^3-(gamm-ph0*h(ix))^3*dhbarphdx(ix)) < 1e-1)
						bslip(ix) =		0.5*(bslip(ix-1)+bslip(ix+1));
					end
				end
			end
		end
	end
		
	p0 =			-(kapp/(gamm^3))*((muvisc*Q)/(L^2))*(barph(nx)*h(nx));
%-----------------------------------------------------------------------------------------------------------------------------------
	
%-----------------------------------------------------------------------------------------------------------------------------------
%	feeding variables to output struct
	inversesoln.bslip =	bslip;
	inversesoln.p0 =	p0;
%-----------------------------------------------------------------------------------------------------------------------------------
end