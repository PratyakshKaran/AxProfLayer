clear;
fclose all;
close all;
clc;

%-----------------------------------------------------------------------------------------------------------------------------------
% setting up the system

% 1st Representative System
%{
system.L1 =			5.0e-2;				% channel left end
system.L2 =			5.0e-2;				% channel right end
system.H =			10.0e-6;			% channel undeformed height
system.Delt =		1.0e-3;				% solid layer thickness
system.Q =			6.25e-10;			% per unit depth volumetric flow rate
system.p0 =			35.0e-2;			% outlet pressure
system.bslipper =	1.0e-2;				% slip length periodicity
%}

% 2nd Representative System
%
system.L1 =			5.0e-1;				% channel left end
system.L2 =			5.0e-1;				% channel right end
system.H =			10.0e-5;			% channel undeformed height
system.Delt =		1.0e-2;				% solid layer thickness
system.Q =			3.7*6.25e-8;		% per unit depth volumetric flow rate
system.p0 =			35.0e-1;			% outlet pressure
system.bslipper =	1.0e-1;				% slip length periodicity
%

% common parameters for both systems
system.muvisc =		1e-3;				% fluid dynamic viscosity
system.rhdens =		1e3;				% fluid density
system.Ey =			9.5e3*((1+0.49)/(1+0.46));				% solid Young's modulus cross-linker variation amplitude
system.nuPois =		0.49;				% solid Poisson's ratio
system.cl_inlet =	15;					% solid inlet cross-linking ratio
system.cl_outlet =	15;					% solid outlet cross-linking ratio
system.bslipamp =	00.0e-6;			% slip length amplitude
system.bslipsteep =	0.20;				% slip length step tanh-smoothening parameter
system.gapmax =		2;					% desired profile for gap

% simulation constants
simul.nx =			1001;				% nodes on x-axis (general, solid domain and fluid domain)
simul.ny =			501;				% nodes on y-axis (fluid domain)
simul.nbary =		501;				% nodes on y-axis (solid domain)
simul.errtol =		1e-6;				% error tolerance for non-linear solver
simul.itermax =		1001;				% maximum number of iterations for non-linear solver
simul.relax =		0.75;				% non-linear solver iteration under-relaxation parameter
simul.flowcalc =	1;					% switch for calculation of fluid domain variables
simul.deformcalc =	1;					% switch for calculation of solid domain variables
simul.casevis =		1;					% switch to generate plots
simul.iguess =		0;					% switch for initial guess value (0 - no guess, 1 - fed as argument, 2 - base solution)
simul.pguess =		0;					% switch for initial guess value (0 - no guess, 1 - fed as argument, 2 - base solution)
simul.npos =		11;					% number of x-positions to show the velocity field at
simul.modeshear =	0;					% which option of shear rate control is wanted (0 - bottom wall, 1 - top wall, 2 - average)
simul.showNRprog =	1;					% show progress of NewtonRaphson solver

% setting the x-axis
system.x =				linspace(-1/(system.bslipper/system.L1),(system.L2/system.L1)/(system.bslipper/system.L1),simul.nx);

% setting elasticity axial variation
Ey =					system.Ey*ones(1,simul.nx);
system.lambd =			((Ey*system.nuPois)/((1+system.nuPois)*(1-2*system.nuPois)));
system.G =				Ey/(2*(1+system.nuPois));

% setting platform axial profile
system.xigreek =		ones(1,simul.nx);

% obtaining parameter values
params =				paramscalc(system,simul);

%{
Inverse Feed In:
Fig 2b:
system.gap =			0.5*(system.gapmax*linspace(2.0,2.0,simul.nx)+1.0* ...
						atanh(((2*(params.alph-params.kapp*0.95*system.x))/(1+params.alph))-1));
trail =					spline(system.x(end-4:end-1),system.gap(end-4:end-1),system.x(end-4:end));
system.bslip =			zeros(1,simul.nx);
Fig 3b:
system.gap =			system.gapmax*linspace(1.0,1.0,simul.nx)-0.95*(system.gapmax-1)* ...
						tanh(((system.bslipper/system.L1)/system.bslipsteep)*(system.x));
system.bslip =			zeros(1,simul.nx);
Fig 4b:
system.gap =			system.gapmax*linspace(1.0,1.0,simul.nx)-0.95*(system.gapmax-1)*cos(pi*linspace(-1.0,1.0,simul.nx));
system.bslip =			zeros(1,simul.nx);
Fig 5b:
system.gap =			system.gapmax*linspace(1.0,1.0,simul.nx)-0.95*(system.gapmax-1)*sin(pi*linspace(-1.0,1.0,simul.nx));
system.bslip =			zeros(1,simul.nx);
Fig 2a:
system.bslip =			system.H*10*(0.5*(1+cos((2*pi*(system.x))/10))).^10;
system.p0 =				218.75;
Fig 3a:
system.bslip =			system.H*10*((0.5*(1+cos((2*pi*(system.x-3.8))/10))).^10+(0.5*(1+cos((2*pi*(system.x+3.8))/10))).^10);
system.p0 =				218.75;
%}

system.gap =			system.gapmax*linspace(1.0,1.0,simul.nx)-0.95*(system.gapmax-1)*cos(pi*linspace(-1.0,1.0,simul.nx));
system.bslip =			zeros(1,simul.nx);

%{ 
for subfigs b in figs 2-5:
system.gap =			<copy paste>;
system.bslip =			zeros(1,simul.nx);
for subfigs a in figs 2-3:
system.bslip =			<copy paste>;
system.p0 =				<copy paste>;
for subfig a in fig 4:
system.bslip =			zeros(1,simul.nx);
system.p0 =				70.0e2;
system.gap =			linspace(1.62,1.5,simul.nx);
for subfig a in fig 5:
system.bslip =			zeros(1,simul.nx);
system.p0 =				70.0e2;
system.gap =			linspace(1.5,1.45,simul.nx);
for subfig a in fig 6:
system.gap =			system.gapmax*linspace(1.0,1.0,simul.nx)-0.95*(system.gapmax-1)* ...
						tanh(((system.bslipper/system.L1)/system.bslipsteep)*(system.x));
system.bslip =			zeros(1,simul.nx);
for subfig b in fig 6:
system.gap =			system.gapmax*linspace(1.0,1.0,simul.nx)-0.95*(system.gapmax-1)*cos(pi*linspace(-1.0,1.0,simul.nx));
system.bslip =			zeros(1,simul.nx);
%}

system.gap =			system.gapmax*linspace(1.0,1.0,simul.nx)-0.95*(system.gapmax-1)*cos(pi*linspace(-1.0,1.0,simul.nx));
system.bslip =			zeros(1,simul.nx);

%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
%{ 
for subfigs b in figs 2-5:
system.h =				(params.gamm/params.ph0)*(1-system.gap);
tapersoln =				tapersolve(system,simul);
system.xigreek =		tapersoln.xigreek;
system.Delt =			tapersoln.Delt;
system.p0 =				tapersoln.p0;
soln =					casesolve(system,simul);
for subfigs a in figs 2-3:
soln =					casesolve(system,simul);
for subfig a in fig 4:
system.h =				(params.gamm/params.ph0)*(1-system.gap);
tapersoln =				tapersolve(system,simul);
system.xigreek =		tapersoln.xigreek;
system.Delt =			tapersoln.Delt;
system.p0 =				tapersoln.p0;
system.bslip =			system.H*10*0.5*(1+tanh(((system.bslipper/system.L1)/system.bslipsteep)*(system.x)));
soln =					casesolve(system,simul);
for subfig a in fig 5:
system.h =				(params.gamm/params.ph0)*(1-system.gap);
tapersoln =				tapersolve(system,simul);
system.xigreek =		tapersoln.xigreek;
system.Delt =			tapersoln.Delt;
system.p0 =				tapersoln.p0;
system.bslip =			system.H*10*((0.5*(1+cos((2*pi*(system.x-3.8))/10))).^10+(0.5*(1+cos((2*pi*(system.x+3.8))/10))).^10);
soln =					casesolve(system,simul);
for subfig a in fig 6:
system.h =				(params.gamm/params.ph0)*(1-system.gap);
system.dvxdy =			0.12*ones(1,simul.nx);
shearsoln =				shearsolve(system,simul);
system.bslip =			shearsoln.bslip;
tapersoln =				tapersolve(system,simul);
system.xigreek =		tapersoln.xigreek;
system.Delt =			tapersoln.Delt;
system.p0 =				tapersoln.p0;
soln =					casesolve(system,simul);
for subfig b in fig 6:
system.h =				(params.gamm/params.ph0)*(1-system.gap);
system.dvxdy =			0.12*ones(1,simul.nx);
shearsoln =				shearsolve(system,simul);
system.bslip =			shearsoln.bslip;
tapersoln =				tapersolve(system,simul);
system.xigreek =		tapersoln.xigreek;
system.Delt =			tapersoln.Delt;
system.p0 =				tapersoln.p0;
soln =					casesolve(system,simul);
%}

system.h =				(params.gamm/params.ph0)*(1-system.gap);
system.dvxdy =			0.12*ones(1,simul.nx);
shearsoln =				shearsolve(system,simul);
system.bslip =			shearsoln.bslip;
tapersoln =				tapersolve(system,simul);
system.xigreek =		tapersoln.xigreek;
system.Delt =			tapersoln.Delt;
system.p0 =				tapersoln.p0;
soln =					casesolve(system,simul);

soln.xigreekspread =	...
						zeros(simul.nx,simul.nbary);
for iy = 1:simul.nbary
	for ix = 1:simul.nx
		system.xigreekspread(ix,iy) =		...
						system.xigreek(ix);
	end
end	
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
% presenting the solution
disp([			'$$ \alpha=',num2str(soln.alph),	',\beta=',num2str(soln.bet),	',\gamma=',num2str(soln.gamm),	...
				',\kappa=',num2str(soln.kapp),',\eta=',num2str(soln.et),		',\phi0=',num2str(soln.ph0),'	$$']);
disp(			['Advection term coefficient is ',num2str(((soln.gamm)/(soln.kapp))*((system.rhdens*system.Q)/system.muvisc))]);

if (simul.casevis ~= 0)
	
	fig1 = figure();
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
	plot(ax1,soln.x,system.bslip*1e6,'r-','LineWidth',3);
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xlabel(ax1,'$x$','Interpreter','Latex','FontSize',25);
	ylabel(ax1,'$b~(\mu$m)','Interpreter','Latex','FontSize',25);
	saveas(fig1,'bslip_x.fig');
	close(fig1);
	
	fig1 = figure();
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
	plot(ax1,soln.x,system.bslip/system.H,'r-','LineWidth',3);
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xlabel(ax1,'$x$','Interpreter','Latex','FontSize',25);
	ylabel(ax1,'$\displaystyle \frac{b}{H}$','Interpreter','Latex','FontSize',25);
	saveas(fig1,'bslipnd_x.fig');
	close(fig1);

	fig1 = figure();
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
	plot(ax1,soln.x,system.xigreek,'r-','LineWidth',3);
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xlabel(ax1,'$x$','Interpreter','Latex','FontSize',25);
	ylabel(ax1,'$\xi$','Interpreter','Latex','FontSize',25);
	saveas(fig1,'xigreek_x.fig');
	close(fig1);

	fig1 = figure();
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
	plot(ax1,soln.x,soln.p,'r-','LineWidth',3);
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xlabel(ax1,'$x$','Interpreter','Latex','FontSize',25);
	ylabel(ax1,'$p$','Interpreter','Latex','FontSize',25);
	saveas(fig1,'p_x.fig');
	close(fig1);

	fig1 = figure();
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
	plot(ax1,soln.x,soln.dpdx(1,:),'r-','LineWidth',3);
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xlabel(ax1,'$x$','Interpreter','Latex','FontSize',25);
	ylabel(ax1,'$\displaystyle \frac{dp}{dx}$','Interpreter','Latex','FontSize',25);
	saveas(fig1,'dpdx_x.fig');
	close(fig1);

	fig1 = figure();
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
	plot(ax1,soln.x,soln.h,'b-','LineWidth',3);
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xlabel(ax1,'$x$','Interpreter','Latex','FontSize',25);
	ylabel(ax1,'$h$','Interpreter','Latex','FontSize',25);
	saveas(fig1,'h_x.fig');
	close(fig1);

	fig1 = figure();
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
	plot(ax1,soln.x,soln.dhdx(1,:),'b-','LineWidth',3);
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xlabel(ax1,'$x$','Interpreter','Latex','FontSize',25);
	ylabel(ax1,'$\displaystyle \frac{dh}{dx}$','Interpreter','Latex','FontSize',25);
	saveas(fig1,'dhdx_x.fig');
	close(fig1);

	
	fig1 = figure();
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
	plot(ax1,soln.x,1-(soln.ph0*soln.h)/soln.gamm,'b-','LineWidth',3);
	hold on;
	plot(ax1,soln.x,1-0*(soln.ph0*soln.h)/soln.gamm,'c:','LineWidth',3);
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xlabel(ax1,'$x$','Interpreter','Latex','FontSize',25);
	ylabel(ax1,'$\displaystyle 1-\frac{\phi h}{\gamma}$ (gap)','Interpreter','Latex','FontSize',25);
	ylims =	get(ax1,'ylim');
	set(ax1,'ylim',[0,ylims(2)]);
	saveas(fig1,'gap_x.fig');
	close(fig1);

	if (simul.casevis == 2)

		fig1 = figure();
		set(fig1,'Position',[33.33,33.33,745,745]);
		ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
		plot(ax1,soln.xd*1e6,soln.pd*1e-6,'r-','LineWidth',3);
		set(ax1,'TickDir','out');
		set(ax1,'XMinorTick','on','YMinorTick','on');
		set(ax1,'FontSize',25);	
		set(ax1,'TickLabelInterpreter','latex');
		set(ax1,'box','on');
		xlabel(ax1,'$x^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
		ylabel(ax1,'$p^{*}$ (MPa)','Interpreter','Latex','FontSize',25);
		saveas(fig1,'pd_xd.fig');
		close(fig1);

		fig1 = figure();
		set(fig1,'Position',[33.33,33.33,745,745]);
		ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
		plot(ax1,soln.xd*1e6,system.bslip*1e6,'r-','LineWidth',3);
		set(ax1,'TickDir','out');
		set(ax1,'XMinorTick','on','YMinorTick','on');
		set(ax1,'FontSize',25);	
		set(ax1,'TickLabelInterpreter','latex');
		set(ax1,'box','on');
		xlabel(ax1,'$x^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
		ylabel(ax1,'$b~(\mu$m)','Interpreter','Latex','FontSize',25);
		saveas(fig1,'bslip_xd.fig');
		close(fig1);

		fig1 = figure();
		set(fig1,'Position',[33.33,33.33,745,745]);
		ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
		plot(ax1,soln.xd*1e6,system.bslip/system.H,'r-','LineWidth',3);
		set(ax1,'TickDir','out');
		set(ax1,'XMinorTick','on','YMinorTick','on');
		set(ax1,'FontSize',25);	
		set(ax1,'TickLabelInterpreter','latex');
		set(ax1,'box','on');
		xlabel(ax1,'$x^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
		ylabel(ax1,'$\displaystyle \frac{b}{H}$','Interpreter','Latex','FontSize',25);
		saveas(fig1,'bslipnd_xd.fig');
		close(fig1);

		fig1 = figure();
		set(fig1,'Position',[33.33,33.33,745,745]);
		ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
		plot(ax1,soln.xd*1e6,system.xigreek,'r-','LineWidth',3);
		set(ax1,'TickDir','out');
		set(ax1,'XMinorTick','on','YMinorTick','on');
		set(ax1,'FontSize',25);	
		set(ax1,'TickLabelInterpreter','latex');
		set(ax1,'box','on');
		xlabel(ax1,'$x^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
		ylabel(ax1,'$\xi$','Interpreter','Latex','FontSize',25);
		saveas(fig1,'xigreek_xd.fig');
		close(fig1);

		fig1 = figure();
		set(fig1,'Position',[33.33,33.33,745,745]);
		ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
		plot(ax1,soln.xd*1e6,soln.dpdxd(1,:),'r-','LineWidth',3);
		set(ax1,'TickDir','out');
		set(ax1,'XMinorTick','on','YMinorTick','on');
		set(ax1,'FontSize',25);	
		set(ax1,'TickLabelInterpreter','latex');
		set(ax1,'box','on');
		xlabel(ax1,'$x^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
		ylabel(ax1,'$\displaystyle \frac{dp^{*}}{dx^*}$ (MPa/$\mu$m)','Interpreter','Latex','FontSize',25);
		saveas(fig1,'dpddxd_xd.fig');
		close(fig1);

		fig1 = figure();
		set(fig1,'Position',[33.33,33.33,745,745]);
		ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
		plot(ax1,soln.xd*1e6,soln.hd*1e6,'b-','LineWidth',3);
		set(ax1,'TickDir','out');
		set(ax1,'XMinorTick','on','YMinorTick','on');
		set(ax1,'FontSize',25);	
		set(ax1,'TickLabelInterpreter','latex');
		set(ax1,'box','on');
		xlabel(ax1,'$x^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
		ylabel(ax1,'$h^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
		saveas(fig1,'hd_xd.fig');
		close(fig1);

		fig1 = figure();
		set(fig1,'Position',[33.33,33.33,745,745]);
		ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
		plot(ax1,soln.xd*1e6,soln.dhdxd(1,:),'b-','LineWidth',3);
		set(ax1,'TickDir','out');
		set(ax1,'XMinorTick','on','YMinorTick','on');
		set(ax1,'FontSize',25);	
		set(ax1,'TickLabelInterpreter','latex');
		set(ax1,'box','on');
		xlabel(ax1,'$x^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
		ylabel(ax1,'$\displaystyle \frac{dh^{*}}{dx^*}$ (m/m)','Interpreter','Latex','FontSize',25);
		saveas(fig1,'dhddxd_xd.fig');
		close(fig1);

		fig1 = figure();
		set(fig1,'Position',[33.33,33.33,745,745]);
		ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
		plot(ax1,soln.xd*1e6,(1-(soln.ph0*soln.h)/soln.gamm)*soln.gamm*soln.L*1e6,'b-','LineWidth',3);
		hold on;
		plot(ax1,soln.xd*1e6,((1-0*(soln.ph0*soln.h)/soln.gamm))*soln.gamm*soln.L*1e6,'c:','LineWidth',3);
		set(ax1,'TickDir','out');
		set(ax1,'XMinorTick','on','YMinorTick','on');
		set(ax1,'FontSize',25);	
		set(ax1,'TickLabelInterpreter','latex');
		set(ax1,'box','on');
		xlabel(ax1,'$x^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
		ylabel(ax1,'$H-h^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
		ylims =	get(ax1,'ylim');
		set(gca,'ylim',[0,ylims(2)]);
		saveas(fig1,'gapd_xd.fig');
		close(fig1);

		if (simul.flowcalc == 1)

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			surface(soln.fluid.xfluid,soln.fluid.yfluid,soln.fluid.dvxdxfluid,'LineStyle','none');
			colorbar();
			hold on;
			plot(ax1,soln.x,1-(soln.ph0*soln.h)/soln.gamm,'b-','LineWidth',3);
			plot(ax1,soln.x,0*(1-(soln.ph0*soln.h)/soln.gamm),'k-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x$','Interpreter','Latex','FontSize',25);
			ylabel('$y$','Interpreter','Latex','FontSize',25);
			title(ax1,'$\displaystyle \frac{dv_x}{dx}$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'extensionalrate.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			plot(soln.fluid.xfluid(1:end,1),soln.fluid.dvxdxfluid(1:end,1),'b-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x$','Interpreter','Latex','FontSize',25);
			ylabel('$\displaystyle \frac{dv_x}{dx}$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'extensionalrate_bottom.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			plot(soln.fluid.xfluid(1:end,end),soln.fluid.dvxdxfluid(1:end,end),'b-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x$','Interpreter','Latex','FontSize',25);
			ylabel('$\displaystyle \frac{dv_x}{dx}$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'extensionalrate_top.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			plot(soln.fluid.xfluid(1:end,(simul.ny-1)/2+1),soln.fluid.dvxdxfluid(1:end,(simul.ny-1)/2+1),'b-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x$','Interpreter','Latex','FontSize',25);
			ylabel('$\displaystyle \frac{dv_x}{dx}$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'extensionalrate_mid.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			surface(soln.fluid.xfluid,soln.fluid.yfluid,soln.fluid.dvxdyfluid,'LineStyle','none');
			colorbar();
			hold on;
			plot(ax1,soln.x,1-(soln.ph0*soln.h)/soln.gamm,'b-','LineWidth',3);
			plot(ax1,soln.x,0*(1-(soln.ph0*soln.h)/soln.gamm),'k-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x$','Interpreter','Latex','FontSize',25);
			ylabel('$y$','Interpreter','Latex','FontSize',25);
			title(ax1,'$\displaystyle \frac{dv_x}{dy}$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'shearrate.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			plot(soln.fluid.xfluid(1:end,1),soln.fluid.dvxdyfluid(1:end,1),'b-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x$','Interpreter','Latex','FontSize',25);
			ylabel('$\displaystyle \frac{dv_x}{dy}$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'shearrate_bottom.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			plot(soln.fluid.xfluid(1:end,end),soln.fluid.dvxdyfluid(1:end,end),'b-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x$','Interpreter','Latex','FontSize',25);
			ylabel('$\displaystyle \frac{dv_x}{dy}$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'shearrate_top.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			plot(soln.fluid.xfluid(1:end,(simul.ny-1)/2+1),soln.fluid.dvxdyfluid(1:end,(simul.ny-1)/2+1),'b-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x$','Interpreter','Latex','FontSize',25);
			ylabel('$\displaystyle \frac{dv_x}{dy}$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'shearrate_mid.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			surface(soln.fluid.xfluidd*1e6,soln.fluid.yfluidd*1e6,	...
					soln.fluid.dvxdxfluidd,'LineStyle','none');
			colorbar();
			hold on;
			plot(ax1,soln.xd*1e6,(1-(soln.ph0*soln.h)/soln.gamm)*soln.gamm*soln.L*1e6,'b-','LineWidth',3);
			plot(ax1,soln.xd*1e6,0*(1-(soln.ph0*soln.h)/soln.gamm)*soln.gamm*soln.L*1e6,'k-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
			ylabel('$y^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
			title(ax1,'$\displaystyle \frac{dv_x^*}{dx^*} (s^{-1})$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'extensionalrated.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			plot(soln.fluid.xfluidd(1:end,1)*1e6,soln.fluid.dvxdxfluidd(1:end,1),'b-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
			ylabel('$\displaystyle \frac{dv_x^*}{dx^*} (s^{-1})$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'extensionalrated_bottom.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			plot(soln.fluid.xfluidd(1:end,end)*1e6,soln.fluid.dvxdxfluidd(1:end,end),'b-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
			ylabel('$\displaystyle \frac{dv_x^*}{dx^*} (s^{-1})$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'extensionalrated_top.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			plot(soln.fluid.xfluidd(1:end,(simul.ny-1)/2+1)*1e6,soln.fluid.dvxdxfluidd(1:end,(simul.ny-1)/2+1),'b-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
			ylabel('$\displaystyle \frac{dv_x^*}{dx^*} (s^{-1})$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'extensionalrated_mid.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			surface(soln.fluid.xfluidd*1e6,soln.fluid.yfluidd*1e6,soln.fluid.dvxdyfluidd,'LineStyle','None');
			colorbar();
			hold on;
			plot(ax1,soln.xd*1e6,(1-(soln.ph0*soln.h)/soln.gamm)*soln.gamm*soln.L*1e6,'b-','LineWidth',3);
			plot(ax1,soln.xd*1e6,0*(1-(soln.ph0*soln.h)/soln.gamm)*soln.gamm*soln.L*1e6,'k-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
			ylabel('$y^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
			title(ax1,'$\displaystyle \frac{dv_x^*}{dy^*} (s^{-1})$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'shearrated.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			plot(soln.fluid.xfluidd(1:end,1)*1e6,soln.fluid.dvxdyfluidd(1:end,1),'b-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
			ylabel('$\displaystyle \frac{dv_x^*}{dy^*} (s^{-1})$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'shearrated_bottom.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			plot(soln.fluid.xfluidd(1:end,end)*1e6,soln.fluid.dvxdyfluidd(1:end,end),'b-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
			ylabel('$\displaystyle \frac{dv_x^*}{dy^*} (s^{-1})$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'shearrated_top.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			plot(soln.fluid.xfluidd(1:end,(simul.ny-1)/2+1)*1e6,soln.fluid.dvxdyfluidd(1:end,(simul.ny-1)/2+1),'b-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x^{*}$ ($\mu$m)','Interpreter','Latex','FontSize',25);
			ylabel('$\displaystyle \frac{dv_x^*}{dy^*} (s^{-1})$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'shearrated_mid.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			for ipos = 1:simul.npos
			plot(ax1,soln.fluid.yfluid((ipos-1)*((simul.nx-1)/(simul.npos-1))+1,1:simul.ny),	...
				soln.fluid.vxfluid((ipos-1)*((simul.nx-1)/(simul.npos-1))+1,1:simul.ny),'LineWidth',3);		
				hold on;
			end
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$y$','Interpreter','Latex','FontSize',25);
			ylabel('$v_x$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'quartervx.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			xin = zeros(1,simul.npos*simul.ny);
			xin1= zeros(1,simul.npos*simul.ny);
			yin = zeros(1,simul.npos*simul.ny);
			cin = zeros(1,simul.npos*simul.ny);
			for ipos = 2:simul.npos-1
				xin((ipos-1)*simul.ny+1:(ipos-1)*simul.ny+simul.ny) = soln.x((ipos-1)*((simul.nx-1)/(simul.npos-1))+1)+ ...
				(soln.fluid.vxfluid((ipos-1)*((simul.nx-1)/(simul.npos-1))+1,1:simul.ny)/max(max(soln.fluid.vxfluid)))* ...
				0.75*(soln.x((3-1)*((simul.nx-1)/(simul.npos-1))+1)-soln.x((2-1)*((simul.nx-1)/(simul.npos-1))+1));
				yin((ipos-1)*simul.ny+1:(ipos-1)*simul.ny+simul.ny) = ...
				soln.fluid.yfluid((ipos-1)*((simul.nx-1)/(simul.npos-1))+1,1:simul.ny);		
				xin1((ipos-1)*simul.ny+1:(ipos-1)*simul.ny+simul.ny) = soln.x((ipos-1)*((simul.nx-1)/(simul.npos-1))+1);
				yin((ipos-1)*simul.ny+1:(ipos-1)*simul.ny+simul.ny) = ...
				soln.fluid.yfluid((ipos-1)*((simul.nx-1)/(simul.npos-1))+1,1:simul.ny);		
				cin((ipos-1)*simul.ny+1:(ipos-1)*simul.ny+simul.ny) = ...
					system.bslip((ipos-1)*((simul.nx-1)/(simul.npos-1))+1)/system.H;
			end
			zin = zeros(size(xin));
			surface([xin;xin],[yin;yin],[zin;zin],[cin;cin],'facecolor','none','edgecolor','interp','linewidth',2,	...
					'LineStyle','none','marker','.');
			colorbar();
			hold on;
			plot(xin1,yin,'.','LineWidth',3,'MarkerSize',1,'color',[0.6,0.6,0.6]);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x$','Interpreter','Latex','FontSize',25);
			ylabel('$y$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'propervx.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			for ipos = 1:simul.npos
			plot(ax1,soln.fluid.yfluidd((ipos-1)*((simul.nx-1)/(simul.npos-1))+1,1:simul.ny),	...
				soln.fluid.vxfluidd((ipos-1)*((simul.nx-1)/(simul.npos-1))+1,1:simul.ny),'LineWidth',3);		
				hold on;
			end
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$y$','Interpreter','Latex','FontSize',25);
			ylabel('$v_x$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'quartervxd.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			xin = zeros(1,simul.npos*simul.ny);
			xin1= zeros(1,simul.npos*simul.ny);
			yin = zeros(1,simul.npos*simul.ny);
			cin = zeros(1,simul.npos*simul.ny);
			for ipos = 2:simul.npos-1
				xin((ipos-1)*simul.ny+1:(ipos-1)*simul.ny+simul.ny) = (soln.xd((ipos-1)*((simul.nx-1)/(simul.npos-1))+1)+ ...
				(soln.fluid.vxfluidd((ipos-1)*((simul.nx-1)/(simul.npos-1))+1,1:simul.ny)/max(max(soln.fluid.vxfluidd)))* ...
				0.75*(soln.xd((3-1)*((simul.nx-1)/(simul.npos-1))+1)-soln.xd((2-1)*((simul.nx-1)/(simul.npos-1))+1)))*1e6;
				yin((ipos-1)*simul.ny+1:(ipos-1)*simul.ny+simul.ny) = ...
				soln.fluid.yfluidd((ipos-1)*((simul.nx-1)/(simul.npos-1))+1,1:simul.ny)*1e6;		
				cin((ipos-1)*simul.ny+1:(ipos-1)*simul.ny+simul.ny) = system.bslip((ipos-1)*((simul.nx-1)/(simul.npos-1))+1);
				xin1((ipos-1)*simul.ny+1:(ipos-1)*simul.ny+simul.ny) = (soln.xd((ipos-1)*((simul.nx-1)/(simul.npos-1))+1))*1e6;
			end
			zin = zeros(size(xin));
			surface([xin;xin],[yin;yin],[zin;zin],[cin;cin],'facecolor','none','edgecolor','interp','linewidth',2,	...
					'LineStyle','none','marker','.');
			colorbar();
			hold on;
			plot(xin1,yin,'.','LineWidth',3,'MarkerSize',1,'color',[0.6,0.6,0.6]);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x^*~(\mu$m)','Interpreter','Latex','FontSize',25);
			ylabel('$y^*~(\mu$m)','Interpreter','Latex','FontSize',25);
			saveas(fig1,'propervxd.fig');
			close(fig1);

		end

		if (simul.deformcalc == 1)

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			surface(soln.solid.xsolid,soln.solid.barysolid,soln.solid.ubarysolid,'LineStyle','none');
			colorbar();
			hold on;
			grid off;
			plot(ax1,soln.x,system.xigreek,'b-','LineWidth',3);
			plot(ax1,soln.x,0*system.xigreek,'k-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x$','Interpreter','Latex','FontSize',25);
			ylabel('$y$','Interpreter','Latex','FontSize',25);
			title(ax1,'$u_{\bar{y}}$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'heatmapdeformLagrangian.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			surface(soln.solid.xsolid,soln.solid.barysolid+(soln.ph0/soln.bet)*soln.solid.ubarysolid,soln.solid.ubarysolid,	...
				'LineStyle','none');
			colorbar();
			hold on;
			grid off;
			plot(ax1,soln.x,system.xigreek,'b-','LineWidth',3);
			plot(ax1,soln.x,-(soln.ph0/soln.bet)*soln.h,'k-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x^*~(\mu$m)','Interpreter','Latex','FontSize',25);
			ylabel('$y^*~(\mu$m)','Interpreter','Latex','FontSize',25);
			title(ax1,'$u_{\bar{y}}^*~(\mu$m)','Interpreter','Latex','FontSize',25);
			saveas(fig1,'heatmapdeformEulerian.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			surface(soln.solid.xsolidd*1e6,soln.solid.barysolidd*1e6,soln.solid.ubarysolidd*1e6,'LineStyle','none');
			colorbar();
			hold on;
			grid off;
			plot(ax1,soln.xd*1e6,system.xigreek*soln.L*soln.bet*1e6,'b-','LineWidth',3);
			plot(ax1,soln.xd*1e6,0*system.xigreek*soln.L*soln.bet*1e6,'k-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x^*~(\mu$m)','Interpreter','Latex','FontSize',25);
			ylabel('$y^*~(\mu$m)','Interpreter','Latex','FontSize',25);
			title(ax1,'$u_{\bar{y}}^*~(\mu$m)','Interpreter','Latex','FontSize',25);
			saveas(fig1,'heatmapdeformLagrangiand.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			surface(soln.solid.xsolidd*1e6,soln.solid.barysolidd*1e6+soln.solid.ubarysolidd*1e6,soln.solid.ubarysolidd*1e6, ...
			'LineStyle','none');
			colorbar();
			hold on;
			grid off;
			plot(ax1,soln.xd*1e6,system.xigreek*soln.L*soln.bet*1e6,'b-','LineWidth',3);
			plot(ax1,soln.xd*1e6,-soln.h*soln.L*soln.ph0*1e6,'k-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x^*~(\mu$m)','Interpreter','Latex','FontSize',25);
			ylabel('$y^*~(\mu$m)','Interpreter','Latex','FontSize',25);
			title(ax1,'$u_{\bar{y}}^*~(\mu$m)','Interpreter','Latex','FontSize',25);
			saveas(fig1,'heatmapdeformEuleriand.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			xin = zeros(1,simul.npos*simul.nx);
			yin = zeros(1,simul.npos*simul.nx);
			yin1= zeros(1,simul.npos*simul.nx);
			cin = zeros(1,simul.npos*simul.nx);
			for ipos = 1:simul.npos
				yin1((ipos-1)*simul.nx+1:(ipos-1)*simul.nx+simul.nx) = ...
				soln.solid.barysolid(1:simul.nx,(ipos-1)*((simul.nbary-1)/(simul.npos-1))+1)';
				yin((ipos-1)*simul.nx+1:(ipos-1)*simul.nx+simul.nx) = ...
				soln.solid.barysolid(1:simul.nx,(ipos-1)*((simul.nbary-1)/(simul.npos-1))+1)'+...
				(soln.solid.ubarysolid(1:simul.nx,(ipos-1)*((simul.nbary-1)/(simul.npos-1))+1))'*(soln.ph0/soln.bet);
				xin((ipos-1)*simul.nx+1:(ipos-1)*simul.nx+simul.nx) = ...
				soln.solid.xsolid(1:simul.nx,(ipos-1)*((simul.nbary-1)/(simul.npos-1))+1);
				cin((ipos-1)*simul.nx+1:(ipos-1)*simul.nx+simul.nx) = system.xigreek(1:simul.nx);
			end
			zin = zeros(size(xin));
			surface([xin;xin],[yin;yin],[zin;zin],[cin;cin],'facecolor','none','edgecolor','interp','linewidth',2,	...
					'LineStyle','none','marker','.');
			colorbar();
			hold on;
			plot(xin,yin1,'.','LineWidth',3,'MarkerSize',1,'color',[0.6,0.6,0.6]);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x$','Interpreter','Latex','FontSize',25);
			ylabel('$y$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'quarterubary.fig');
			close(fig1);

			fig1 = figure();
			set(fig1,'Position',[33.33,33.33,745,745]);
			ax1 =   axes(fig1,'Position',[0.300,0.300,0.600,0.600]);
			xin = zeros(1,simul.npos*simul.nx);
			yin = zeros(1,simul.npos*simul.nx);
			yin1= zeros(1,simul.npos*simul.nx);
			cin = zeros(1,simul.npos*simul.nx);
			for ipos = 1:simul.npos
				yin1((ipos-1)*simul.nx+1:(ipos-1)*simul.nx+simul.nx) = ...
				soln.solid.barysolidd(1:simul.nx,(ipos-1)*((simul.nbary-1)/(simul.npos-1))+1)';
				yin((ipos-1)*simul.nx+1:(ipos-1)*simul.nx+simul.nx) = ...
				soln.solid.barysolidd(1:simul.nx,(ipos-1)*((simul.nbary-1)/(simul.npos-1))+1)'+...
				(soln.solid.ubarysolidd(1:simul.nx,(ipos-1)*((simul.nbary-1)/(simul.npos-1))+1))';
				xin((ipos-1)*simul.nx+1:(ipos-1)*simul.nx+simul.nx) = ...
				soln.solid.xsolid(1:simul.nx,(ipos-1)*((simul.nbary-1)/(simul.npos-1))+1);
				cin((ipos-1)*simul.nx+1:(ipos-1)*simul.nx+simul.nx) = system.xigreek(1:simul.nx);
			end
			zin = zeros(size(xin));
			surface([xin;xin],[yin;yin],[zin;zin],[cin;cin],'facecolor','none','edgecolor','interp','linewidth',2,	...
					'LineStyle','none','marker','.');
			colorbar();
			hold on;
			plot(xin,yin1,'.','LineWidth',3,'MarkerSize',1,'color',[0.6,0.6,0.6]);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel('$x$','Interpreter','Latex','FontSize',25);
			ylabel('$y$','Interpreter','Latex','FontSize',25);
			saveas(fig1,'quarterubary.fig');
			close(fig1);

		end

	end
		
	if ((simul.deformcalc == 1)	&& (simul.flowcalc == 1))
		
		set(groot,'DefaultFigureColormap',hot(64));
		fig1 = figure();
		set(fig1,'Position',[100,35,750,950]);		
		if (sum(abs(system.bslip)) ~= 0)
			ax0 =   axes(fig1,'Position',[0.225,0.125,0.700,0.175]);
			ax1 =   axes(fig1,'Position',[0.225,0.300,0.700,0.175]);		
			ax2 =   axes(fig1,'Position',[0.225,0.475,0.700,0.250]);		
			ax2a =  axes(fig1,'Position',[0.225,0.475,0.700,0.250]);		
			ax3 =   axes(fig1,'Position',[0.225,0.725,0.700,0.250]);		
			plot(ax0,soln.x,system.bslip/system.H,'b-','LineWidth',3);
			set(ax0,'TickDir','out');
			set(ax0,'XMinorTick','on','YMinorTick','on');
			set(ax0,'FontSize',25);	
			set(ax0,'TickLabelInterpreter','latex');
			set(ax0,'box','on');
			xlabel(ax0,'$x$','Interpreter','Latex','FontSize',25);
			ylabel(ax0,'$\displaystyle \frac{b}{H}$','Interpreter','Latex','FontSize',25);
			set(ax0,'ylim',[min(system.bslip/system.H)-0.15*abs(min(system.bslip/system.H))-1.5,	...
				max(system.bslip/system.H)+0.05*abs(max(system.bslip/system.H))]);
			xlim0 = get(ax0,'xlim');
			plot(ax1,soln.x,1-(soln.ph0*soln.h)/soln.gamm,'k-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			set(ax1,'xtick',[]);
			ylabel(ax1,'$\displaystyle 1-\frac{\phi h}{\gamma}$','Interpreter','Latex','FontSize',25);
			ylims =	get(ax1,'ylim');
			set(ax1,'ylim',[min(1-(soln.ph0*soln.h)/soln.gamm)-0.005*abs(min(1-(soln.ph0*soln.h)/soln.gamm)),	...
				max(1-(soln.ph0*soln.h)/soln.gamm)+0.005*abs(max(1-(soln.ph0*soln.h)/soln.gamm))]);
			set(ax1,'xlim',xlim0);
			xin = zeros(1,simul.npos*simul.ny);
			xin1= zeros(1,simul.npos*simul.ny);
			yin = zeros(1,simul.npos*simul.ny);
			for ipos = 1:simul.npos
				xin((ipos-1)*simul.ny+1:(ipos-1)*simul.ny+simul.ny) = soln.x((ipos-1)*((simul.nx-1)/(simul.npos-1))+1)+ ...
				(soln.fluid.vxfluid((ipos-1)*((simul.nx-1)/(simul.npos-1))+1,1:simul.ny)/max(max(soln.fluid.vxfluid)))* ...
				0.75*(soln.x((3-1)*((simul.nx-1)/(simul.npos-1))+1)-soln.x((2-1)*((simul.nx-1)/(simul.npos-1))+1));
				yin((ipos-1)*simul.ny+1:(ipos-1)*simul.ny+simul.ny) = ...
				soln.fluid.yfluid((ipos-1)*((simul.nx-1)/(simul.npos-1))+1,1:simul.ny);		
				xin1((ipos-1)*simul.ny+1:(ipos-1)*simul.ny+simul.ny) = soln.x((ipos-1)*((simul.nx-1)/(simul.npos-1))+1);
				yin((ipos-1)*simul.ny+1:(ipos-1)*simul.ny+simul.ny) = ...
				soln.fluid.yfluid((ipos-1)*((simul.nx-1)/(simul.npos-1))+1,1:simul.ny);		
			end
			surface(ax2,soln.fluid.xfluid,soln.fluid.yfluid,soln.fluid.dvxdyfluid,'LineStyle','none');
			clbr = colorbar('peer',ax2,'north');
			set(clbr,'ticklabelinterpreter','latex');
			set(ax2,'TickDir','out');
			set(ax2,'XMinorTick','on','YMinorTick','on');
			set(ax2,'FontSize',25);	
			set(ax2,'TickLabelInterpreter','latex');
			set(ax2,'box','on');
			set(ax2,'ycolor','none');
			set(ax2,'xtick',[]);
			set(ax2,'ytick',[]);
			set(ax2,'ylim',[0,max(1-(soln.ph0*soln.h)/soln.gamm)*1.5]);
			set(ax2,'xlim',xlim0);		
			set(fig1,'CurrentAxes',ax2a);
			plot(ax2a,xin1,yin,'.','LineWidth',2,'MarkerSize',1,'color',[0.3,0.3,0.3]);
			hold on;
			plot(ax2a,xin,yin,'.','LineWidth',2,'MarkerSize',6,'color','c');
			set(ax2a,'xtick',[]);
			set(ax2a,'color','none'); set(ax2a,'xcolor','none');
			set(ax2a,'ylim',[0,max(1-(soln.ph0*soln.h)/soln.gamm)*1.5]);
			set(ax2a,'xlim',xlim0);
			set(ax2a,'FontSize',25);	
			ylabel(ax2a,'$y$','Interpreter','Latex','FontSize',25);	
			set(ax2a,'TickLabelInterpreter','latex');
			srf = surface(ax3,soln.solid.xsolid,soln.solid.barysolid,	...
				(soln.ph0/soln.bet)*(soln.solid.ubarysolid./system.xigreekspread),'LineStyle','none');
			hold on;
			grid off;
			set(ax3,'ylim',[0,1.35]);
			clbr = colorbar('peer',ax3,'north');
			set(clbr,'ticklabelinterpreter','latex');
			set(ax3,'TickDir','out');
			set(ax3,'XMinorTick','on','YMinorTick','on');
			set(ax3,'FontSize',25);	
			set(ax3,'TickLabelInterpreter','latex');
			set(ax3,'box','on');
			ylabel(ax3,'$\bar{y}$','Interpreter','Latex','FontSize',25);
			set(ax3,'xlim',xlim0);
			set(ax3,'xtick',[]);
		else
			ax1 =   axes(fig1,'Position',[0.225,0.125,0.700,0.270]);		
			ax2 =   axes(fig1,'Position',[0.225,0.395,0.700,0.290]);		
			ax2a =  axes(fig1,'Position',[0.225,0.395,0.700,0.290]);		
			ax3 =   axes(fig1,'Position',[0.225,0.685,0.700,0.290]);		
			plot(ax1,soln.x,1-(soln.ph0*soln.h)/soln.gamm,'k-','LineWidth',3);
			set(ax1,'TickDir','out');
			set(ax1,'XMinorTick','on','YMinorTick','on');
			set(ax1,'FontSize',25);	
			set(ax1,'TickLabelInterpreter','latex');
			set(ax1,'box','on');
			xlabel(ax1,'$x$','Interpreter','Latex','FontSize',25);
			ylabel(ax1,'$\displaystyle 1-\frac{\phi h}{\gamma}$','Interpreter','Latex','FontSize',25);
			ylims =	get(ax1,'ylim');
			set(ax1,'ylim',[min(1-(soln.ph0*soln.h)/soln.gamm)-0.005*abs(min(1-(soln.ph0*soln.h)/soln.gamm)),	...
				max(1-(soln.ph0*soln.h)/soln.gamm)+0.005*abs(max(1-(soln.ph0*soln.h)/soln.gamm))]);
			xlim0 = get(ax1,'xlim');
			xin = zeros(1,simul.npos*simul.ny);
			xin1= zeros(1,simul.npos*simul.ny);
			yin = zeros(1,simul.npos*simul.ny);
			for ipos = 1:simul.npos
				xin((ipos-1)*simul.ny+1:(ipos-1)*simul.ny+simul.ny) = soln.x((ipos-1)*((simul.nx-1)/(simul.npos-1))+1)+ ...
				(soln.fluid.vxfluid((ipos-1)*((simul.nx-1)/(simul.npos-1))+1,1:simul.ny)/max(max(soln.fluid.vxfluid)))* ...
				0.75*(soln.x((3-1)*((simul.nx-1)/(simul.npos-1))+1)-soln.x((2-1)*((simul.nx-1)/(simul.npos-1))+1));
				yin((ipos-1)*simul.ny+1:(ipos-1)*simul.ny+simul.ny) = ...
				soln.fluid.yfluid((ipos-1)*((simul.nx-1)/(simul.npos-1))+1,1:simul.ny);		
				xin1((ipos-1)*simul.ny+1:(ipos-1)*simul.ny+simul.ny) = soln.x((ipos-1)*((simul.nx-1)/(simul.npos-1))+1);
				yin((ipos-1)*simul.ny+1:(ipos-1)*simul.ny+simul.ny) = ...
				soln.fluid.yfluid((ipos-1)*((simul.nx-1)/(simul.npos-1))+1,1:simul.ny);		
			end
			surface(ax2,soln.fluid.xfluid,soln.fluid.yfluid,soln.fluid.dvxdyfluid,'LineStyle','none');
			clbr = colorbar('peer',ax2,'north');
			set(clbr,'ticklabelinterpreter','latex');
			set(ax2,'TickDir','out');
			set(ax2,'XMinorTick','on','YMinorTick','on');
			set(ax2,'FontSize',25);	
			set(ax2,'TickLabelInterpreter','latex');
			set(ax2,'box','on');
			set(ax2,'ycolor','none');
			set(ax2,'xtick',[]);
			set(ax2,'ytick',[]);
			set(ax2,'ylim',[0,max(1-(soln.ph0*soln.h)/soln.gamm)*1.5]);
			set(ax2,'xlim',xlim0);		
			set(fig1,'CurrentAxes',ax2a);
			plot(ax2a,xin1,yin,'.','LineWidth',2,'MarkerSize',1,'color',[0.3,0.3,0.3]);
			hold on;
			plot(ax2a,xin,yin,'.','LineWidth',2,'MarkerSize',6,'color','c');
			set(ax2a,'xtick',[]);
			set(ax2a,'color','none'); set(ax2a,'xcolor','none');
			set(ax2a,'ylim',[0,max(1-(soln.ph0*soln.h)/soln.gamm)*1.5]);
			set(ax2a,'xlim',xlim0);
			set(ax2a,'FontSize',25);	
			ylabel(ax2a,'$y$','Interpreter','Latex','FontSize',25);	
			set(ax2a,'TickLabelInterpreter','latex');
			srf = surface(ax3,soln.solid.xsolid,soln.solid.barysolid,	...
					(soln.ph0/soln.bet)*(soln.solid.ubarysolid./system.xigreekspread),'LineStyle','none');
			hold on;
			grid off;
			set(ax3,'ylim',[0,1.35]);
			clbr = colorbar('peer',ax3,'north');
			set(clbr,'ticklabelinterpreter','latex');
			set(ax3,'TickDir','out');
			set(ax3,'XMinorTick','on','YMinorTick','on');
			set(ax3,'FontSize',25);	
			set(ax3,'TickLabelInterpreter','latex');
			set(ax3,'box','on');
			ylabel(ax3,'$\bar{y}$','Interpreter','Latex','FontSize',25);
			set(ax3,'xlim',xlim0);
			set(ax3,'xtick',[]);
		end
		saveas(fig1,'system_4panel.fig');
		close(fig1);
	end
	
end
%-----------------------------------------------------------------------------------------------------------------------------------
% OLD:
%{
% setting slip profile
system.bslip =			system.bslipamp*0.5*(1+tanh(((system.bslipper/system.L1)/system.bslipsteep)*(system.x)));
0 system.gap =			system.gapmax*linspace(1.0,1.0,simul.nx);
1 system.gap =			system.gapmax*linspace(1.0,1.0,simul.nx)+0.95*(system.gapmax-1)*sin(pi*linspace(-1.0,1.0,simul.nx));
2 system.gap =			system.gapmax*linspace(1.0,1.0,simul.nx)-0.95*(system.gapmax-1)*sin(pi*linspace(-1.0,1.0,simul.nx));
system.bslip =			...
system.H*10*((0.5*(1+cos((2*pi*(system.x-3.8))/10))).^10+(0.5*(1+cos((2*pi*(system.x+3.8))/10))).^10);
3 system.gap =			system.gapmax*linspace(1.0,1.0,simul.nx)+0.95*(system.gapmax-1)*cos(pi*linspace(-1.0,1.0,simul.nx));
4 system.gap =			system.gapmax*linspace(1.0,1.0,simul.nx)-0.95*(system.gapmax-1)*cos(pi*linspace(-1.0,1.0,simul.nx));
system.bslip =			system.bslipamp*0.5*(1+tanh(((system.bslipper/system.L1)/system.bslipsteep)*(system.x)));
system.dvxdy =			pchip(soln.fluid.xfluid(:,1),soln.fluid.dvxdyfluid(:,1),system.x)*0 + ...
						min(pchip(soln.fluid.xfluid(:,1),soln.fluid.dvxdyfluid(:,1),system.x))*ones(1,simul.nx)*1;
5 system.gap =			system.gapmax*linspace(1.0,1.0,simul.nx)+0.95*(system.gapmax-1)* ...
						tanh(((system.bslipper/system.L1)/system.bslipsteep)*(system.x));
6 system.gap =			system.gapmax*linspace(1.0,1.0,simul.nx)-0.95*(system.gapmax-1)* ...
						tanh(((system.bslipper/system.L1)/system.bslipsteep)*(system.x));
system.bslip =			...
system.H*10*((0.5*(1+cos((2*pi*(system.x-3.8))/10))).^10+(0.5*(1+cos((2*pi*(system.x+3.8))/10))).^10);
7 system.gap =			0.5*(system.gapmax*linspace(2.0,2.0,simul.nx)+1.0* ...
						atanh(((2*(params.alph-params.kapp*0.95*system.x))/(1+params.alph))-1));
trail =					spline(system.x(end-4:end-1),system.gap(end-4:end-1),system.x(end-4:end));
system.gap(end) =		trail(end);
system.bslip =			system.H*10*(0.5*(1+cos((2*pi*(system.x))/10))).^10;
8 system.gap =			0.5*(system.gapmax*linspace(2.0,2.0,simul.nx)+1.0* ...
						atanh(-(((2*(params.alph-params.kapp*0.95*system.x))/(1+params.alph))-1)));
trail =					spline(system.x(end-4:end-1),system.gap(end-4:end-1),system.x(end-4:end));
system.gap(end) =		trail(end);
%}

% solution based on shear rate
%{
system.dvxdy =			pchip(soln.fluid.xfluid(:,1),soln.fluid.dvxdyfluid(:,1),system.x)*0 + ...
						min(pchip(soln.fluid.xfluid(:,1),soln.fluid.dvxdyfluid(:,1),system.x))*ones(1,simul.nx)*1.0;
shearsoln =				shearsolve(system,simul);
system.bslip =			shearsoln.bslip;
tapersoln =				tapersolve(system,simul);
system.xigreek =		tapersoln.xigreek;
system.Delt =			tapersoln.Delt;
system.p0 =				tapersoln.p0;
soln =					casesolve(system,simul);
%}

%{
9 (botched) system.gap =			system.gapmax*linspace(1,1,simul.nx)-0.95*(system.gapmax-1)* ...
						tanh(((system.bslipper/system.L1)/system.bslipsteep)*(system.x));
%}