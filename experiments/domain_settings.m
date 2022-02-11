
% Domain settings

% --------------------------- Spatial dimensions --------------------------

% Rectangular domain
dim.z.min =    0.0;      % Depth increases in the positive z-direction,
dim.z.max =  500.0;      % so z_min = 0 is at the surface of the domain
                                
dim.x.min =    0.0;
dim.x.max = 1000.0;

% ----------------------------- Spatial grid ------------------------------

% Number of grid points along each spatial direction
N.z =  51;
N.x = 101;

% -------------------------- Boundary conditions --------------------------

% 'D' Dirichlet, 'N' Neumann, 'A' Absorbing (i.e. unbounded) 

                                  BCs.z.min = 'A';
%                           --------------------------- 
           BCs.x.min = 'A';                             BCs.x.max = 'A';
%                           --------------------------- 
                                  BCs.z.max = 'A';

% ---------------- Absorbing boundary layer (ABL) settings ----------------

% Thicknesses and strengths of ABLs in terms of number of grid nodes

                              BCs.z.ABL.min.thickness = 10;
                              BCs.z.ABL.min.strength  = 200;
%                             ------------------------------ 
BCs.x.ABL.min.thickness = 10;                               BCs.x.ABL.max.thickness = 10;
BCs.x.ABL.min.strength  = 200;                              BCs.x.ABL.max.strength  = 200;
%                             ------------------------------
                              BCs.z.ABL.max.thickness = 10;
                              BCs.z.ABL.max.strength  = 200;
                            
% -------------------------------------------------------------------------
