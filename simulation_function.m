function [resJ_end,resA_end,end_pop,strain_totalJ,strain_totalA,indexJ_end,indexA_end,RESJ,RESA,DISPREV,NVEC] = simulation_function(t_max,a0,g0,c1a,c2a,c1g,c2g,beta0,alpha,resJmin,resJmax,resJ_start,resAmin,resAmax,resA_start,h,f,resJmutationprob,init_pop,strain_totalJ,strain_totalA,indexJ_start,indexA_start,res0,nevol,version)

% This function runs an evolutionary simulation for a fixed number of
% timesteps.

eqtol = 1e-3;
exttol = 1e-5;

% Set up vectors to be used later:
ResJ = linspace(resJmin,resJmax,res0);
ResA = linspace(resAmin,resAmax,res0);
RESJ = zeros(nevol,res0);
RESA = zeros(nevol,res0);
DISPREV = zeros(nevol,1);
NVEC = zeros(nevol,1);

% Initial conditions:
resJ_current = resJ_start;
resA_current = resA_start;
indexJ_current = indexJ_start;
indexA_current = indexA_start;

% Each value of ievol is one evolutionary timestep.
for ievol=1:nevol
    
    % Find the ecological equilibrium:
    if version==1
        [~,SJ1,SA1,IJ1,IA1,~] = eco_dynamics_function_v1(t_max,a0,g0,c1a,c2a,c1g,c2g,beta0,alpha,resJ_current,resA_current,h,f,eqtol,init_pop,strain_totalJ,strain_totalA);
    elseif version==2
        [~,SJ1,SA1,IJ1,IA1,~] = eco_dynamics_function_v2(t_max,a0,g0,c1a,c2a,c1g,c2g,beta0,alpha,resJ_current,resA_current,h,f,eqtol,init_pop,strain_totalJ,strain_totalA);
    elseif version==3
        [~,SJ1,SA1,IJ1,IA1,~] = eco_dynamics_function_v3(t_max,a0,g0,c1a,c2a,c1g,c2g,beta0,alpha,resJ_current,resA_current,h,f,eqtol,init_pop,strain_totalJ,strain_totalA);
    elseif version==4
        [~,SJ1,SA1,IJ1,IA1,~] = eco_dynamics_function_v4(t_max,a0,g0,c1a,c2a,c1g,c2g,beta0,alpha,resJ_current,resA_current,h,f,eqtol,init_pop,strain_totalJ,strain_totalA);
    elseif version==5
        [~,SJ1,SA1,IJ1,IA1,~] = eco_dynamics_function_v5(t_max,a0,g0,c1a,c2a,c1g,c2g,beta0,alpha,resJ_current,resA_current,h,f,eqtol,init_pop,strain_totalJ,strain_totalA);
    elseif version==6
        [~,SJ1,SA1,IJ1,IA1,~] = eco_dynamics_function_v6(t_max,a0,g0,c1a,c2a,c1g,c2g,beta0,alpha,resJ_current,resA_current,h,f,eqtol,init_pop,strain_totalJ,strain_totalA);
    end
        
    % Re-format this output into a single matrix:
    SJ=zeros(strain_totalJ,strain_totalA);
    SA=zeros(strain_totalJ,strain_totalA);
    IJ=zeros(strain_totalJ,strain_totalA);
    IA=zeros(strain_totalJ,strain_totalA);
    for j=1:strain_totalJ
        for k=1:strain_totalA
            SJ(j,k) = SJ1(end,j+(k-1)*strain_totalJ);
            SA(j,k) = SA1(end,j+(k-1)*strain_totalJ);
            IJ(j,k) = IJ1(end,j+(k-1)*strain_totalJ);
            IA(j,k) = IA1(end,j+(k-1)*strain_totalJ);
        end
    end
    N = SA+SJ+IA+IJ;
    
    % Remove extinct classes
    Nrows=sum(N,2);
    Ntotal=sum(N,'all');
    
    % See if any resJ strains go extinct:
    extinct = (Nrows/Ntotal)<exttol;
    strain_totalJ = strain_totalJ-sum(extinct);
    SA(extinct,:) = [];
    SJ(extinct,:) = [];
    IA(extinct,:) = [];
    IJ(extinct,:) = [];
    N(extinct,:) = [];
    indexJ_current(extinct) = [];
    resJ_current(extinct) = [];
    
    % Update N:
    Ncolumns=sum(N,1);
    Ntotal=sum(N,'all');
    
    % See if any resA strains go extinct:
    extinct1 = (Ncolumns/Ntotal)<exttol;
    strain_totalA = strain_totalA-sum(extinct1);
    SA(:,extinct1) = [];
    SJ(:,extinct1) = [];
    IA(:,extinct1) = [];
    IJ(:,extinct1) = [];
    N(:,extinct1) = [];
    indexA_current(extinct1) = [];
    resA_current(extinct1) = [];
    
    % Update tracker
    Nrows=sum(N,2);
    Ncolumns=sum(N,1);
    Ntotal=sum(N,'all'); 
    
    % Proportion of individuals of each resJ strain
    RESJ(ievol,indexJ_current) = Nrows./Ntotal;
    % Proportion of individuals of each resA strain
    RESA(ievol,indexA_current) = Ncolumns./Ntotal;
    % Proportion of individuals who have the disease
    DISPREV(ievol) = (sum(IA,'all')+sum(IJ,'all'))/Ntotal;
    % Total population density:
    NVEC(ievol) = Ntotal;
    
    Nvector=zeros(4*strain_totalJ*strain_totalA,1);
    for j=1:strain_totalJ
        for k=1:strain_totalA
            Nvector(j+(k-1)*strain_totalJ)=N(j,k);
        end
    end
    
    % If the mutation occurs in the juvenile resistance:
    if (rand<resJmutationprob)
        weightedprob = Nvector/sum(Nvector);
        cumsum1 = cumsum(weightedprob);
        r1 = rand*cumsum1(end);
        mutator_loc = (find(r1<cumsum1,1));
        mutator_locJ=0;
        for j=1:strain_totalJ
            for k=1:strain_totalA
                if mutator_loc==j+(k-1)*strain_totalJ
                   mutator_locJ=j;
                end
            end
        end
        mutator = indexJ_current(mutator_locJ);
 
        if(mutator==1) % Mutate up
            mutant = mutator+1;
        elseif(mutator==res0) % Mutate down
            mutant = mutator-1;
        else
            if(rand>0.5) % Mutate up
                mutant = mutator+1;
            else % Mutate down
                mutant = mutator-1;
            end
        end
        if(~ismember(mutant,indexJ_current)) % New strain
            strain_totalJ = strain_totalJ+1;
            % Update vector of resJ trait values:
            resJ_current_again=NaN(1,length(resJ_current)+1);
            for i=1:length(resJ_current)
                resJ_current_again(1,i)=resJ_current(i);
            end
            resJ_current_again(1,end)=ResJ(mutant);
            resJ_current=resJ_current_again;
            % Update vector of resJ trait value indices:
            indexJ_current_again=NaN(1,length(indexJ_current)+1);
            for i=1:length(indexJ_current)
                indexJ_current_again(1,i)=indexJ_current(i);
            end
            indexJ_current_again(1,end)=mutant;
            indexJ_current=indexJ_current_again;
            % Add a small population with the new trait value:
            SJ_again=NaN(size(SJ,1)+1,size(SJ,2));
            for i=1:size(SJ,1)
                for j=1:size(SJ,2)
                    SJ_again(i,j)=SJ(i,j);
                end
            end
            SJ_again(end,:)=SJ(mutator_locJ,:)/10;
            SJ=SJ_again;
            SA_again=NaN(size(SA,1)+1,size(SA,2));
            for i=1:size(SA,1)
                for j=1:size(SA,2)
                    SA_again(i,j)=SA(i,j);
                end
            end
            SA_again(end,:)=SA(mutator_locJ,:)/10;
            SA=SA_again;
            IJ_again=NaN(size(IJ,1)+1,size(IJ,2));
            for i=1:size(IJ,1)
                for j=1:size(IJ,2)
                    IJ_again(i,j)=IJ(i,j);
                end
            end
            IJ_again(end,:)=IJ(mutator_locJ,:)/10;
            IJ=IJ_again;
            IA_again=NaN(size(IA,1)+1,size(IA,2));
            for i=1:size(IA,1)
                for j=1:size(IA,2)
                    IA_again(i,j)=IA(i,j);
                end
            end
            IA_again(end,:)=IA(mutator_locJ,:)/10;
            IA=IA_again;
            
        end
    end
    
    % If the mutation occurs in the adult resistance:
    if (rand>resJmutationprob) 
        weightedprob = Nvector/sum(Nvector);
        cumsum1 = cumsum(weightedprob);
        r1 = rand*cumsum1(end);
        mutator_loc = (find(r1<cumsum1,1));
        mutator_locA=0;
        for j=1:strain_totalJ
            for k=1:strain_totalA
                if mutator_loc==j+(k-1)*strain_totalJ
                   mutator_locA=k;
                end
            end
        end
        mutator = indexA_current(mutator_locA);
    
        if(mutator==1) % Mutate up
            mutant = mutator+1;
        elseif(mutator==res0) % Mutate down
            mutant = mutator-1;
        else
            if(rand>0.5) % Mutate up
                mutant = mutator+1;
            else % Mutate down
                mutant = mutator-1;
            end
        end
        if(~ismember(mutant,indexA_current)) % New strain
            strain_totalA = strain_totalA+1;
            % Update vector of resA trait values:
            resA_current_again=NaN(1,length(resA_current)+1);
            for i=1:length(resA_current)
                resA_current_again(1,i)=resA_current(i);
            end
            resA_current_again(1,end)=ResA(mutant);
            resA_current=resA_current_again;
            % Update vectors of resA trait value indices:
            indexA_current_again=NaN(1,length(indexA_current)+1);
            for i=1:length(indexA_current)
                indexA_current_again(1,i)=indexA_current(i);
            end
            indexA_current_again(1,end)=mutant;
            indexA_current=indexA_current_again;
            % Add small populations with the new trait value to the
            % population:
            SJ_again=NaN(size(SJ,1),size(SJ,2)+1);
            for i=1:size(SJ,1)
                for j=1:size(SJ,2)
                    SJ_again(i,j)=SJ(i,j);
                end
            end
            SJ_again(:,end)=SJ(:,mutator_locA)/10;
            SJ=SJ_again;
            SA_again=NaN(size(SA,1),size(SA,2)+1);
            for i=1:size(SA,1)
                for j=1:size(SA,2)
                    SA_again(i,j)=SA(i,j);
                end
            end
            SA_again(:,end)=SA(:,mutator_locA)/10;
            SA=SA_again;
            IJ_again=NaN(size(IJ,1),size(IJ,2)+1);
            for i=1:size(IJ,1)
                for j=1:size(IJ,2)
                    IJ_again(i,j)=IJ(i,j);
                end
            end
            IJ_again(:,end)=IJ(:,mutator_locA)/10;
            IJ=IJ_again;
            IA_again=NaN(size(IA,1),size(IA,2)+1);
            for i=1:size(IA,1)
                for j=1:size(IA,2)
                    IA_again(i,j)=IA(i,j);
                end
            end
            IA_again(:,end)=IA(:,mutator_locA)/10;
            IA=IA_again;
            
        end
    end
    
    % Update initial conditions
    init_pop = NaN(1,4*strain_totalJ*strain_totalA);
    for i=1:strain_totalJ*strain_totalA
        init_pop(4*i)=IA(i);
        init_pop(4*i-1)=IJ(i);
        init_pop(4*i-2)=SA(i);
        init_pop(4*i-3)=SJ(i);
    end
end

% Create outputs:
resJ_end=resJ_current;
resA_end=resA_current;
end_pop=init_pop;
indexJ_end=indexJ_current;
indexA_end=indexA_current;

end