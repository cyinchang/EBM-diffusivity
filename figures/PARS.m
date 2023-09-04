function [par] = PARS(parstr, options);
% PARS  Parameters for plots and analyses.

  %%%%%%%%%%%%%           parse options         %%%%%%%%%%%%%%%%%%%%
  error(nargchk(1, 2, nargin))     
  if nargin ==1 | isempty(options)
    fopts      = [];
  else
    fopts      = fieldnames(options);
  end
  
  if strmatch('slide', fopts)
    slide = options.slide;
  else 
    slide = 0;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
  % flag for file dump of plots, and directory for file dumps
  plot_file           = 0;
  plot_dir            = '/home/tapio/tex/dry_scaling/figs/'; 
    
  % plot parameters
  if (slide)
    % size of figures (in inches)
    figure_width      = 4;
    figure_height     = figure_width / 1.6;
    
    % line widths
%    line_width        = 1.0;
    line_width        = 0.5*2;
    thin_line_width   = 0.5;
    axes_line_width   = 1.0;
    thick_line_width  = 1.0;

    % line color
    pos_line_color    = 'r';
    neg_line_color    = 'c';
    line_color        = 'c';
    
    % line styles
    line_style        = '-';
    thin_line_cstyle  = 'm-';
    thick_line_cstyle = 'r-';

    % tick length
    tick_length       = [0.01 0.01];
    
    % marker size
    marker_size       = 10;
    
    % font sizes
    font_size         = 10;
    clabel_size       = 10;
  else
    % size of figures (in inches)

    % Science (one column)
    %figure_width      = 2.3;
    %figure_height     = figure_width;
    
    % Princeton UP
    %figure_width      = 4.5;  % PUP page
    %figure_width      = 2.5;
    %figure_width      = 3.5;
    %figure_height     = figure_width / 1.666;
    
    % AMS journals
    figure_width      = 3.125;
    figure_height     = figure_width / 1.618;

    % Small square 
    %figure_width      = 2.5;
    %figure_height     = figure_width;
    
    % supercriticality figure in SW06 (JAS)
    %figure_width      = 3;
    %figure_height     = figure_width;
    
    % line colors and styles	
    pos_line_color    = 'k';
    neg_line_color    = 'm';
    line_color        = 'k';
    line_style        = '-';
    thin_line_cstyle  = 'k:';
    thick_line_cstyle = 'k-';

    % Annual Reviews
    %figure_width       = 6.5;
    %figure_width      = 3.38;
    %figure_height     = figure_width/1.6;
    %pos_line_color    = 'm-';
    %neg_line_color    = 'm--';
    %line_color        = 'm';
    %thin_line_cstyle  = 'c';
    %thick_line_cstyle = 'r-';

    % line widths
    line_width        = .75;  % pt 
    thin_line_width   = .5; % PUP
    axes_line_width   = 1;
    thick_line_width  = 1.2;

    % tick length
    tick_length       = [0.025 0.025];
    
    % marker size
    marker_size       = 4.5;
    
    % font sizes
    font_size         = 7;
    clabel_size       = 7;
  end
  
  % physical parameters
  radius_earth 	      = 6.3710e6;     % in meters
  radius 	      = 6.3710e6;     % in meters
  gravity 	      = 9.80;         % in m/s^2
  gas_constant	      = 287.04;       % gas constant dry air [J/kg/K] 
  kappa               = 2/7;
  cp                  = gas_constant/kappa; % specific heat at constant pressure 
  cp_v                = 1870;         % specific heat water vapor [J/kg/K]
  cl                  = 4190;         % heat capacity liquid water [J/kg/K]
  gas_constant_v      = 461.50;       % gas constant water vapor [J/kg/K]
  omega               = 7.292e-5;     % angular velocity earth [rad/s]
  mean_sfc_press      = 1e5;          % mean surface pressure
  l_cond              = 2.5e6;        % latent heat cond [J/kg] at zero deg C
  latent_heat_v       = 2.5e6;        % latent heat of evaporation [J/kg]
  stefan              = 5.6734e-8;    % Stefan-Boltzmann constant [W/m^2/K^4]


  gray                = [.5 .5 .5];   % RGB value for gray lines

  % for conversion of units
  secs_per_day        = 86400;        % length of Earth day [seconds]
  deg                 = pi/180;
  
  % for isentropic-coordinate plots
  sfc_cdf_isolines    = [.5 .5];      % for plotting of median only
  pressure_surfaces  = [250 500 750 925]; %isobars to be plotted
  
  eval(['par = ', parstr, ';']);
  
  
