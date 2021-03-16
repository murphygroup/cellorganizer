function [result, success] = generate_walk_path(model, shape_space, energies,param)
% [result, success] = generate_walk_path(...
%   y2, convex_hull, tes, param.save_filename)
%
% Generate a random walk path through a shape space (not the shapes, points in the shape space corresponding to them) specified by y2, convex_hull, and tes using various random walk models ('brownian', 'constantstep', 'willmore', 'density').
%
%Edited:
%6/26/13 D. Sullivan - added directed walk using willmore energy and
%                      created switch statement instead of if-elseif
%                      for efficiency.

%error('Implementation yet unfinished below this line!')
if isempty(energies)
    if strcmpi(param.walk_type,'willmore')
        warning(['No energies specified, assuming uniform energies ',...
            'and therefor defaulting to "brownian" walk.']);
        param.walk_type = 'brownian';
    end
end

[y2, convex_hull, tes] = shape_space{:};
result = struct();
result.walk_path = [];
result.start_index = [];

number_images = size(y2, 1);
numdims = size(y2, 2);
walk_path = zeros(0, numdims);
start_index = 0;

success = false;
% Check if work has been done
[can_start, final_name, final_exists] = ...
    chunk_start(param.save_filename);
if (~can_start)
    stack = dbstack();
    warning('Skipping "%s" at line %d in file %s.', ...
        param.save_filename, stack(1).line, stack(1).file )
    if final_exists
        load([param.save_filename '.mat'])
        
        % If there is a saved random state, set the current state to
        % it for identical downstream random number generation:
        if (exist('current_random_state', 'var'))
            RandStream.setGlobalStream(current_random_state)
        end
        
        success = true;
    else
        %error('Cannot continue, serial 1 being computed elsewhere.')
        error('  Being computed elsewhere.')
    end
else
    if ~isfield(param,'start_index')
        start_index = ceil(rand(1) .* number_images)
    else
        start_index = param.start_index;
    end
    %find the neighborhood
    D = pdist2(y2,y2);
    [distance,index] = sort(D,2);
    current_point = y2(start_index, :);
    %D. Sullivan 6/26/13 added energy tracking for Willmore mode
    current_energy = energies(start_index,:);
    energypred = zeros(1,param.mmax);
    predicted_energies = zeros(1,param.mmax);
    %     if strcmpi(param.walk_type,'density')
    %         N = size(y2,1);
    %         d = size(y2,2);
    %         H = BW_Estimation(y2');
    %         Hd = diag(H);
    %     end
    if strcmpi(param.walk_type,'lowess') || strcmpi(param.walk_type,'loess')
        param.maxiter = param.mmax;
        %Currently this method is only supported in 2D fitting
        current_point = current_point(:,1:2);
        minvals = min(y2);
        maxvals = max(y2);
        param.range = [minvals(1),maxvals(1),minvals(2),maxvals(2)];
        [xopt,fopt,niter,gnorm,walk_path,fvalue] = graddirectedwalk(current_point,model.dynamic,param)
    else
        
        for m = 1:param.mmax
            walk_path(m,:)=current_point;
            
            while( true )
                %D. Sullivan 6/26/13 added switch statement instead of if-elseif
                %to allow for easier and cleaner expansion to other walk types
                switch lower(param.walk_type)
                    case {'brownian'}
                        candidate = current_point + randn(1, numdims) .* ...
                            param.Dc * sqrt(param.dt);
                    case {'constantstep'}
                        direction = zeros( 1, numdims );
                        magnitude = 0;
                        while( magnitude == 0 )
                            direction = randn( 1, numdims );
                            magnitude = sqrt(sum( direction.^2 ));
                        end
                        
                        direction = direction./magnitude;
                        candidate = current_point+direction.*param.Dc*param.dt;
                    case {'willmore','density'}
                        %randomly pick a candidate point
                        candidate = current_point + randn(1, numdims) .* ...
                            param.Dc * sqrt(param.dt);
                        %estimate the canditate's Willmore energy using linear
                        %interpolation between neighboring points
                        
                        
                        %first find the nearest neighbors
                        [idx,d] = knnsearch(y2,candidate,'K',3);
                        %                     [idx,d] = knnsearch(energies,current_energy,'K',3);
                        %%%%%%%%%%
                        D =pdist2(y2(idx,:),y2(idx,:));
                        interpE1 = energies(idx(1))+(energies(idx(2))...
                            -energies(idx(1)))*(d(1)/D(1,2));
                        interpE2 = energies(idx(1))+(energies(idx(3))...
                            -energies(idx(1)))*(d(1)/D(1,3));
                        interpE3 = energies(idx(2))+(energies(idx(3))...
                            -energies(idx(2)))*(d(2)/D(2,3));
                        candidateE = mean([interpE1,interpE2,interpE3]);
                        %%%%%%%%%%
                        %                     candidateE = interp3(y2(idx,1),y2(idx,2),y2(idx,3),...
                        %                         energies(idx),candidate(1),candidate(2),...
                        %                         candidate(3));
                        
                        %                     candidateE = abs((energies(idx(1))-energies(idx(2)))/d(2)...
                        %                         +(energies(idx(1))-energies(idx(3)))/d(3));
                        %
                        %now choose whether to accept the move based on energy
                        alpha = current_energy/candidateE;
                        if alpha>=1
                            %The move is energetically favorable, always accept
                            current_point = candidate;
                            current_energy = candidateE;
                            
                            %once we have accepted a move, add it to the energy
                            %space to make it findable in the knn search
                            %%%D. Sullivan 6/29/13 - need to add check if the
                            %%%point is already in the space or we will bias
                            %%%the interpolation/energy guess.
                            if find(energies==current_energy)
                                energies(end+1) = current_energy;
                                y2(end+1,:) = current_point;
                            end
                            
                        else
                            %The move costs energy, sample acceptance
                            %D. Sullivan 6/26/13 - here the probability
                            %distribution of acceptance could be learned from
                            %real data, but for now default to linear
                            %probability.
                            %                         fraction = alpha/(1+exp(-alpha));
                            Paccept = alpha;
                            acceptCrit = rand;%uniformly random select a number 0-1
                            %If my fraction is larger than my acceptCrit,
                            %acccept the move
                            if Paccept>acceptCrit
                                
                                %                             current_point = candidate;
                                current_energy = candidateE;
                                
                                %once we have accepted a move, add it to the
                                %energy space to make it findable in the knn
                                %search
                                %%%D. Sullivan 6/29/13 - need to add check if the
                                %%%point is already in the space or we will bias
                                %%%the interpolation/energy guess.
                                if find(energies==current_energy)
                                    energies(end+1) = current_energy;
                                    y2(end+1,:) = current_point;
                                end
                                
                            else
                                candidate = current_point;
                                %                             current_energy = current_energy;
                                %In Metropolis-Hastings we stay here, not sure
                                %that is realistic in our case. May want to
                                %resample until a move is accepted.
                                %                             current_point = current_point;
                                %                             current_energy = current_energy;
                                %no need to add energy, because we did not
                                %move.
                            end
                            
                        end
                        
                        predicted_energies(m) = current_energy;
                        
                    otherwise
                        error('Unrecognized walk type.');
                end
                if inhull( candidate, y2, convex_hull )
                    current_point = candidate;
                    break;
                    %             elseif strcmpi(param.walk_type,'willmore')||strcmpi(param.walk_type,'density')
                    %                 current_point = current_point;
                    %                 break;
                end
            end
            
        end
    end
    %D. Sullivan 6/26/13 changed to getGlobalStream based on MATLAB warning
    current_random_state = RandStream.getGlobalStream;
    %     current_random_state = RandStream.getDefaultStream;
    save([param.save_filename '.mat'], 'walk_path', 'start_index', 'current_random_state');
    %D. Sullivan 6/29/13 added save energy predictions
    if strcmpi(param.walk_type,'willmore')||strcmpi(param.walk_type,'density')
        save(['temp' filesep 'predicted_energies.mat'],'predicted_energies')
    end
    %%%%%%%%%%%%%%%%
    
    % Save work
    chunk_finish('.', param.save_filename);
end

result.walk_path =walk_path ;
result.start_index = start_index;

success = true;
