This ReadMe file contains information on files associated with the paper: "The evolution of age-specific resistance to infectious disease" by Buckingham, Bruns & Ashby.

All files are provided "as is", without any express or implied warranty. The authors accept no liability for damages of any kind. 

Author: Lydia Buckingham, University of Bath
Contact: ljb74@bath.ac.uk
Date: 06/09/22

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
Figure 1:   "tradeoff_plots.m"
Figure 2:   "f.fig.m"
Figure 3:   "phaseplane_fig_mixedparams.m"
Figure 4:   "beta0_fig_misedparams.m"
Figure 5:   "alpha_fig.m"
Figure S1:  "phaseplane_fig.m"
Figure S2:  "beta0_fig.m"
Figure S3:  "beta0_fig.m"                    - f = 0.3
Figure S4:  "alpha_fig_3curves.m"
Figure S5:  "alpha_fig_3curves.m"            - a0 = 5   &   g0vec    = [0.5,1,5]
Figure S6:  "alpha_fig_3curves.m"            - a0 = 5   &   beta0vec = [5,8,20]
Figure S7:  "alpha_fig_3curves.m"            - a0 = 5   &   c1vec    = [0.3,0.5,0.7]
Figure S8:  "alpha_fig_3curves.m"            - a0 = 5   &   c2vec    = [-2,-3,-8]
Figure S9:  "f_fig_3curves.m"
Figure S10: "f_fig_3curves.m"                - a0 = 5   &   g0vec    = [0.5,1,5]
Figure S11: "f_fig_3curves.m"                - a0 = 5   &   beta0vec = [5,8,20]
Figure S12: "f_fig_3curves.m"                - a0 = 5   &   c1vec    = [0.3,0.5,0.7]
Figure S13: "f_fig_3curves.m"                - a0 = 5   &   c2vec    = [-2,-3,-8]

Note that "alpha_fig_3curves" and "f_fig_3curves" each plot three curves corresponding to three different values of a particular parameter. This varying parameter is different for the different 
supplementary figures and so 'a0vec' should be replaced with 'a0' and the varying trait with its corresponding vector to create these figures. Where 'a0' is assigned a value from 'a0vec' at three 
places in the code, this will also need to be changed to assign the new varying parameter to a value from its vector. 

In all code, 'c2g' is used instead of '-c2J' and 'c2a' is used instead of '-c2A'. 

In all code, the parameter 'h' represents the relative transmissibility of disease from juveniles. This parameter was cut from the model and so is always equal to one in the code.

MEX files written in C# must be compiled before use, using " mex codename.c ".

Be aware that some of the code relies on numerical approximations and so results may not be exact. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DESCRIPTION OF FILES

Figure 1:
"tradeoff_plots.m"		 - Source code for plotting example trade-off functions										 - written in matlab (R2019b). 



Figure 3:
"phaseplane_fig_mixedparams.m"	 - Source code for plotting multiple phase planes using different parameter values 						 - written in matlab (R2019b). 

"sim_inout_function.m" 		 - Function for determining the endpoint of an evolutionary trajectory given initial conditions					 - written in matlab (R2019b). 
"sim_traj_function.m" 		 - Function for simulating an evolutionary trajecotory given initial conditions							 - written in matlab (R2019b). 
"simulation_function.m"          - Function which runs an evolutionary simulation for a fixed number of timesteps						 - written in matlab (R2019b). 
"eco_dynamics_function_v1.c" 	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using 'version 1' trade-offs 	 - written as a MEX file in C#. 
"eco_dynamics_function_v2.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using 'version 2' trade-offs 	 - written as a MEX file in C#. 
"eco_dynamics_function_v3.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using 'version 3' trade-offs  	 - written as a MEX file in C#. 
"eco_dynamics_function_v4.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using 'version 4' trade-offs 	 - written as a MEX file in C#. 
"eco_dynamics_function_v5.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using 'version 5' trade-offs 	 - written as a MEX file in C#. 
"eco_dynamics_function_v6.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using 'version 6' trade-offs 	 - written as a MEX file in C#. 
"fitness_gradients.m"		 - Function for calculating the two fitness gradients at different values of juvenile and adult resistance			 - written in matlab (R2019b). 
"endemic_equilibrium_function.m" - Function for determining the endemic equilibrium of the ecological system							 - written in matlab (R2019b). 
"fast_ode_function.c"            - Function for running a system of ODEs to find its equilibrium								 - written as a MEX file in C#. 



Figure S1:
"phaseplane_fig.m"               - Source code for plotting phase planes                                                                                         - written in matlab (R2019b). 

"fitness_gradients.m"		 - Function for calculating the two fitness gradients at different values of juvenile and adult resistance			 - written in matlab (R2019b). 
"endemic_equilibrium_function.m" - Function for determining the endemic equilibrium of the ecological system							 - written in matlab (R2019b). 
"fast_ode_function.c"            - Function for running a system of ODEs to find its equilibrium								 - written as a MEX file in C#. 


All other figures:
"beta0_fig.m"	                 - Source code for plotting the effect of baseline transmissibility on singular strategies of juvenile and adult resistance	 - written in matlab (R2019b). 
"f_fig.m"		         - Source code for plotting the effect of sterility virulence on singular strategies of juvenile and adult resistance		 - written in matlab (R2019b). 
"alpha_fig.m"	                 - Source code for plotting the effect of mortality virulence on singular strategies of juvenile and adult resistance  		 - written in matlab (R2019b). 
"beta0_fig_mixedparams.m"        - Source code for plotting the effect of baseline transmissibility on singular strategies for different parameter values	 - written in matlab (R2019b). 
"f_fig_3curves.m"		 - Source code for plotting the effect of sterility virulence on singular strategies for 3 different values of another parameter - written in matlab (R2019b). 
"alpha_fig_3curves.m"	         - Source code for plotting the effect of mortality virulence on singular strategies for 3 different values of another parameter - written in matlab (R2019b). 

"beta0_fig_data.m"		 - Function for calculating singular strategies for different values of baseline transmissibility				 - written in matlab (R2019b). 
"f_fig_data.m"			 - Function for calculating singular strategies for different values of sterility virulence					 - written in matlab (R2019b). 
"alpha_fig_data.m"		 - Function for calculating singular strategies for different values of mortality virulence					 - written in matlab (R2019b). 
"fitness_gradients.m"		 - Function for calculating the two fitness gradients at different values of juvenile and adult resistance			 - written in matlab (R2019b). 
"fitgrad_signchange_function.m"	 - Function which determines when the fitness gradient changes sign							  	 - written in matlab (R2019b). 
"singstrats_at_0or1.m"		 - Function which determines singular strategies where one of the resistance traits takes a value of zero or one		 - written in matlab (R2019b). 
"singstrat_finder_analytical.m"	 - Function for calculating singular strategies analytically									 - written in matlab (R2019b). 
"classification_function.m"      - Function for determining whether a singular strategy is an evolutionary repeller, a CSS or a branching point                  - written in matlab (R2019b). 
"endemic_equilibrium_function.m" - Function for determining the endemic equilibrium of the ecological system							 - written in matlab (R2019b). 
"fitgrad_functions.m"		 - Function which determines analytical expressions for fitness gradients and their derivatives					 - written in matlab (R2019b). 
"ode_fast.c"			 - Function for running a system of ODEs to find its equilibrium								 - written as a MEX file in C#. 



Other code included:
"ecological_stability.m"         - Code which shows that there is a unique, linearly stable ecological equilibrium (when R0>1) for a wide range of parameters    - written in matlab (R2019b). 
"simulation_in_steps.m"          - Code which plots an evolutionary trajectory of juvenile and adult resistance                                                  - written in matlab (R2019b). 


See code for full description and instructions for use. 