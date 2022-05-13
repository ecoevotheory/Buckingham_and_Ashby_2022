This ReadMe file contains information on files associated with the paper: "The evolution of age-specific resistance to infectious disease" by Buckingham & Ashby.

All files are provided "as is", without any express or implied warranty. The authors accept no liability for damages of any kind. 

Author: Lydia Buckingham, University of Bath
Contact: ljb74@bath.ac.uk
Date: 11/05/22

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COMMENTS

All code refers to different combinations of trade-offs as numbered 'versions' of the model, as follows:
version = 1 : Juvenile resistance trades off with maturation 	     - Adult resistance trades off with reproduction
version = 2 : Juvenile resistance trades off with juvenile mortality - Adult resistance trades off with adult mortality
version = 3 : Juvenile resistance trades off with reproduction 	     - Adult resistance trades off with adult mortality
version = 4 : Juvenile resistance trades off with maturation	     - Adult resistance trades off with adult mortality
version = 5 : Juvenile resistance trades off with juvenile mortality - Adult resistance trades off with reproduction
version = 6 : Juvenile resistance trades off with reproduction	     - Adult resistance trades off with reproduction

Figures are generated using the source code below and parameter values as given in the code with the following exceptions:
Figure 1: "tradeoff_plots.m"
Figure 2: "effect_of_f_figure.m" & "f_fig_data.m"
Figure 3: "effect_of_beta0_figure.m" & "beta0_fig_data.m"  - a0 = 3
Figure 4: "effect_of_beta0_figure.m" & "beta0_fig_data.m"  -  f = 0.3
Figure 5: "effect_of_alpha_figure.m" & "alpha_fig_data.m"
Figure S1: "effect_of_f_figure.m" & "f_fig_data.m"	   - g0 = 3
Figure S2: "effect_of_beta0_figure.m" & "beta0_fig_data.m" - a0 = 3.5
Figure S3: "effect_of_beta0_figure.m" & "beta0_fig_data.m" - g0 = 0.5
Figure S4: "effect_of_alpha_figure.m" & "alpha_fig_data.m" - g0 = 5
Figure S5: "phaseplane_figure.m"

In all code, 'c2g' is used instead of '-c2J' and 'c2a' is used instead of '-c2A'. 

In all code, the parameter 'h' represents the relative transmissibility of disease from juveniles. This parameter was cut from the model and so is always equal to one in the code.

MEX files written in C# must be compiled before use, using " mex codename.c ".

Be aware that much of the code relies on numerical approximations and so results may not be exact. 
In particular, disease prevalence may not be calculated precisely when the total population density is very low.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DESCRIPTION OF FILES

Figure 1:
"tradeoff_plots.m"		 - Source code for plotting example trade-off functions										 - written in matlab (R2019b). 


Figure S5:
"phaseplane_figure.m"		 - Source code for plotting a phase plane											 - written in matlab (R2019b). 

"sim_inout_function.m" 		 - Function for determining the endpoint of an evolutionary trajectory given initial conditions					 - written in matlab (R2019b). 
"sim_traj_function.m" 		 - Function for simulating an evolutionary trajecotory given initial conditions							 - written in matlab (R2019b). 
"simulation_in_steps_function.m" - Function which runs an evolutionary simulation for a fixed number of timesteps						 - written in matlab (R2019b). 
"plantmodel_function_mex.c" 	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using 'version 1' trade-offs 	 - written as a MEX file in C#. 
"plantmodel_function_v2_mex.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using 'version 2' trade-offs 	 - written as a MEX file in C#. 
"plantmodel_function_v3_mex.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using 'version 3' trade-offs  	 - written as a MEX file in C#. 
"plantmodel_function_v4_mex.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using 'version 4' trade-offs 	 - written as a MEX file in C#. 
"plantmodel_function_v5_mex.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using 'version 5' trade-offs 	 - written as a MEX file in C#. 
"plantmodel_function_v6_mex.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using 'version 6' trade-offs 	 - written as a MEX file in C#. 
"fitness_gradients.m"		 - Function for calculating the two fitness gradients at different values of juvenile and adult resistance			 - written in matlab (R2019b). 
"endemic_equilibrium_function.m" - Function for determining the endemic equilibrium of the ecological system							 - written in matlab (R2019b). 
"ode_fast.c"			 - Function for running a system of ODEs to find its equilibrium								 - written as a MEX file in C#. 


All other figures:
"effect_of_f_figure.m"		 - Source code for plotting the effect of sterility virulence on singular strategies of juvenile and adult resistance		 - written in matlab (R2019b). 
"effect_of_beta0_figure.m"	 - Source code for plotting the effect of baseline transmissibility on singular strategies of juvenile and adult resistance	 - written in matlab (R2019b). 
"effect_of_alpha_figure.m"	 - Source code for plotting the effect of mortality virulence on singular strategies of juvenile and adult resistance  		 - written in matlab (R2019b). 

"f_fig_data.m"			 - Function for calculating singular strategies for different values of sterility virulence					 - written in matlab (R2019b). 
"beta0_fig_data.m"		 - Function for calculating singular strategies for different values of baseline transmissibility				 - written in matlab (R2019b). 
"alpha_fig_data.m"		 - Function for calculating singular strategies for different values of mortality virulence					 - written in matlab (R2019b). 
"fitness_gradients.m"		 - Function for calculating the two fitness gradients at different values of juvenile and adult resistance			 - written in matlab (R2019b). 
"fitgrad_signchange_function.m"	 - Function which determines when the fitness gradient changes sign							  	 - written in matlab (R2019b). 
"singstrats_at_0or1.m"		 - Function which determines singular strategies where one of the resistance traits takes a value of zero or one		 - written in matlab (R2019b). 
"singstrat_finder_analytical.m"	 - Function for calculating singular strategies analytically									 - written in matlab (R2019b). 
"endemic_equilibrium_function.m" - Function for determining the endemic equilibrium of the ecological system							 - written in matlab (R2019b). 
"fitgrad_functions.m"		 - Function which determines analytical expressions for fitness gradients and their derivatives					 - written in matlab (R2019b). 
"ode_fast.c"			 - Function for running a system of ODEs to find its equilibrium								 - written as a MEX file in C#. 


See code for full description and instructions for use. 