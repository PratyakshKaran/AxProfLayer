clear;
fclose all;
close all;
clc;

%-----------------------------------------------------------------------------------------------------------------------------------
% setting up the system
%-----------------------------------------------------------------------------------------------------------------------------------
system.L1 =			5.0e-2;				% channel left end
system.L2 =			5.0e-2;				% channel right end
system.H =			10.0e-6;			% channel undeformed height
system.Delt =		1.0e-3;				% solid layer thickness
system.Q =			1.0e-4;				% per unit depth volumetric flow rate
system.p0 =			35.0e2;				% outlet pressure
system.muvisc =		1e-3;				% fluid dynamic viscosity
system.rhdens =		1e3;				% fluid density
system.Ey =			9.5e3;				% solid Young's modulus cross-linker variation amplitude
system.nuPois =		0.46;				% solid Poisson's ratio
system.cl_inlet =	15;					% solid inlet cross-linking ratio
system.cl_outlet =	15;					% solid outlet cross-linking ratio
system.bslipamp =	00.0e-6;			% slip length amplitude
system.bslipper =	1.0e-2;				% slip length periodicity
system.bslipsteep =	0.25;				% slip length step tanh-smoothening parameter
system.gapmax =		2;					% desired profile for gap

simul.nx =			1001;				% nodes on x-axis (general, solid domain and fluid domain)
simul.ny =			501;				% nodes on y-axis (fluid domain)
simul.nbary =		501;				% nodes on y-axis (solid domain)
simul.errtol =		1e-8;				% error tolerance for non-linear solver
simul.itermax =		1001;				% maximum number of iterations for non-linear solver
simul.relax =		0.75;				% non-linear solver iteration under-relaxation parameter
simul.flowcalc =	1;					% switch for calculation of fluid domain variables
simul.deformcalc =	1;					% switch for calculation of solid domain variables
simul.casevis =		1;					% switch to generate plots
simul.iguess =		0;					% switch for initial guess value (0 - no guess, 1 - fed as argument, 2 - base solution)
simul.pguess =		0;					% switch for initial guess value (0 - no guess, 1 - fed as argument, 2 - base solution)
simul.npos =		11;					% number of x-positions to show the velocity field at
simul.modeshear =	0;					% which option of shear rate control is wanted (0 - bottom wall, 1 - top wall, 2 - average)
simul.nstrong =		21;				% points of strength parametric sweep

%

% setting the x-axis
system.x =				linspace(-1/(system.bslipper/system.L1),(system.L2/system.L1)/(system.bslipper/system.L1),simul.nx);

% setting elasticity axial variation
Ey =					system.Ey*ones(1,simul.nx);
system.lambd =			((Ey*system.nuPois)/((1+system.nuPois)*(1-2*system.nuPois)));
system.G =				Ey/(2*(1+system.nuPois));

% setting platform axial profile
% system.xigreek =		0.5+linspace(0.6,1.0,simul.nx)-(0.5*(1+cos(2*pi*linspace(0,1,simul.nx)-pi)));
% setting slip profile
system.bslip =			zeros(1,simul.nx);
params =				params(system,simul);
system.gap =			system.gapmax*linspace(1.0,1.0,simul.nx)-0.95*(system.gapmax-1)*cos(pi*linspace(-1.0,1.0,simul.nx));
system.h =				(params.gamm/params.ph0)*(1-system.gap);
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------------------
% obtaining the taper inverse solution
tapersoln =				tapersolve(system,simul);
system.Delt =			tapersoln.Delt;
system.p0 =				tapersoln.p0;
% obtaining the sweep
strong =			logspace(-2,0,simul.nstrong);
xparam =			zeros(simul.nstrong,simul.nx);
taperparam =		zeros(simul.nstrong,simul.nx);
gapparam =			zeros(simul.nstrong,simul.nx);
for istrong = 1:simul.nstrong
	system.xigreek =		ones(1,simul.nx)+ ...
							strong(istrong)*(tapersoln.xigreek-ones(1,simul.nx));
	soln =					casesolve(system,simul);
%	presenting the solution
	disp([			'$$ \alpha=',num2str(soln.alph),	',\beta=',num2str(soln.bet),	',\gamma=',num2str(soln.gamm),	...
					',\kappa=',num2str(soln.kapp),',\eta=',num2str(soln.et),		',\phi0=',num2str(soln.ph0),'	$$']);
	disp(			['Advection term coefficient is ',num2str(((soln.gamm)/(soln.kapp))*((system.rhdens*system.Q)/system.muvisc))]);
	for ix = 1:simul.nx
		gapparam(istrong,ix) =		1-(soln.ph0/soln.gamm)*soln.h(ix);
		xparam(istrong,ix) =		soln.x(ix);
		taperparam(istrong,ix) =	strong(istrong);
	end
end
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------

if (simul.casevis == 1)

	fig1 = figure();
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
	wt = waterfall(ax1,xparam,taperparam,gapparam);
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on','ZMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	set(ax1,'yscale','log');
	xlabel(ax1,'$x$','Interpreter','Latex','FontSize',25);
	ylabel(ax1,'taper','Interpreter','Latex','FontSize',25);
	zlabel(ax1,'$\displaystyle 1-\frac{\phi_0 h}{\gamma}$','Interpreter','Latex','FontSize',25);
	saveas(fig1,'param_taper.fig');
	zlim = get(ax1,'zlim');
	set(ax1,'zlim',[0,zlim(2)*1.05]);
	hold on;	
	close(fig1);

end
%-----------------------------------------------------------------------------------------------------------------------------------
