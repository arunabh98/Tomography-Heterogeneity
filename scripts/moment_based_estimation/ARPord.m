function [thetasestimated] = ...
    ARPord(Pgiven, svector, sigmaNoise, initialThetaEstimate,...
           noisyOrientations, angle_amplitude)

    addpath(genpath('../utilities'));

    numkeep = size(Pgiven, 2);
    kmax = numkeep - 1;
    if noisyOrientations == 1
        numstarts = 2;
    else
        numstarts = 12;
    end
    Ord = 8;

    refinedProjections = Pgiven;

    numiter = 40;   

    thetasestimated_bystart = zeros(numkeep, numstarts);
    Ehlccvalues_bystart = zeros(numstarts, 1);

    if noisyOrientations == 1
        min_limit = initialThetaEstimate - angle_amplitude;
        max_limit = initialThetaEstimate + angle_amplitude;
    else
        min_limit = repmat(-179, numkeep, 1);
        max_limit = repmat(180, numkeep, 1);
    end

    parfor start = 1:numstarts

        PMord = ...
            assemblePMord(Pgiven, kmax, svector, Ord, numkeep);
        
        if noisyOrientations == 1
            thetasestimated = initialThetaEstimate;
        else
            % Initialize angles randomly
            thetasestimated = randi([1 179], numkeep, 1);
        end

        A = assembleA(thetasestimated, Ord); 
        IMestimated = A \ PMord;

        Ehlccvec = A * IMestimated - PMord; 
        Ehlcc = norm(Ehlccvec, 2);
        Ehlccvalues = zeros(numiter * (numkeep - 1), 1);
        deltas = zeros(numiter, 1);
        ctr = 1;

        Ehlccvalues(ctr) = Ehlcc;

        for iteration = 1:numiter
            disp([num2str(start), ', ', num2str(iteration)]);

            for i = 1:numkeep
                ctr = ctr + 1;

                besttheta = thetasestimated(i);
                bestE = Ehlcc;

                max_limit_i = max_limit(i);
                min_limit_i = min_limit(i);
                 
                err_t = err_for_all_angles(...
                    Pgiven, kmax, svector, Ord, thetasestimated,...
                    max_limit_i, min_limit_i, i, PMord);
                
                [E_t, idx_err_t] = min(err_t);
                if E_t < bestE
                    bestE = E_t;
                    besttheta = idx_err_t + min_limit(i) - 1;
                end        

                thetasestimated(i) = besttheta;

                A = assembleA(thetasestimated, Ord);
                IMestimated = A \ PMord;

                Ehlccvec = A * IMestimated - PMord;
                Ehlcc = norm(Ehlccvec, 2);
                Ehlccvalues(ctr) = Ehlcc;

            end
            Ehlcciter = Ehlccvalues(ctr);
            Ehlccpreviter = Ehlccvalues(ctr - (numkeep)); % -1 iff theta0 is not updated

            delta = Ehlccpreviter - Ehlcciter;
            deltas(iteration) = delta;
            if abs(delta) < 0.01
                break
            end

        end
        Ehlccvalues_bystart(start) = Ehlccvalues(ctr);
        thetasestimated_bystart(:, start) = thetasestimated;
    end % multistart loop

    [~, optstart] = min(Ehlccvalues_bystart);
    thetasestimated = thetasestimated_bystart(:, optstart);
end
