function TAG_DMOEA(Global)
% <algorithm> <A>
% A twin-archive guided decomposition based multi/many-objective evolutionary algorithm

%------------------------------- Reference --------------------------------
% Sri Srinivasa Raju M, Rammohan Mallipeddi, Kedar Nath Das,A twin-archive
% guided decomposition based multi/many-objective evolutionary algorithm, 
% Swarm and Evolutionary Computation, Volume 71, 2022, 101082, ISSN 2210-6502, 
% https://doi.org/10.1016/j.swevo.2022.101082.

%--------------------------------------------------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate the reference points and random population
    [W,Global.N]        = UniformPoint(Global.N,Global.M);
    Population          = Global.Initialization();
    [z,znad]            = deal(min(Population.objs),max(Population.objs));
    EI                  = zeros(1,Global.N);
    Archive1            = [];
    Archive2            = [];
    wa                  =  0;
    ZNad(Global.gen,:)  = znad;
     %% Neighbourhood
    T               = ceil(Global.N/10);
    B               = pdist2(W,W);
    [~,B]           = sort(B,2);
    B               = B(:,1:T);
    %% Optimization
    while Global.NotTermination(Population) 
    %% Mating Selection and Offspring Generation
        MatingPool                        = randi(Global.N,1,Global.N);  
        Offspring                         = GA(Population(MatingPool));
    %% Environmental Selection tDEA
       [Population,z,znad,EI]             = EnvironmentalSelectiontDEA([Population,Offspring],W,Global.N,z,znad,EI,(Global.gen/Global.maxgen)>=0.5);
    %% Archive Selection       
       [Archive1,Archive2,W,EI,znad,B]    = ArchiveSelectionbased2REA([Archive1,Population],[Archive2,Population,Offspring],W,Global.N,EI,znad,B,wa,~mod(Global.gen,ceil(0.04*Global.maxgen)),z);
                 ZNad(Global.gen,:)       = znad;
                      wa                  = WeightAdapationIndicator(wa,Global,ZNad);
    end
end