function [meas, x_true, phiEst, raylMean, A_sp, nFx, nFy, pupDom] = precompute(dataname)
    % Inputs:
    %   dataname:       name of dataset

    dataSet = load(dataname);
    meas = dataSet.meas;
    x_true = dataSet.z_coef;
    phiEst = dataSet.phiEst;
    raylMean = dataSet.raylMean;
    if nargout == 4
        return
    end
    %% Define optimization parameters
    nFx             = dataSet.nFx;
    nFy             = dataSet.nFy;
    U               = dataSet.U;
    V               = dataSet.V;
    pupDom          = dataSet.pupDom;
    n_img           = length(U);
    MsCell          = cell(n_img,1);
    nonzeroDesign   = cell(n_img,1);
    %% Create Zernike basis and design matrix
    maxDeg  = 5;
    [B,Z]   = zernikeBasis( maxDeg, nFx, nFy ); 
    BMat    = reshape(B,[],size(B,3));
    
    %% Compute the (sparse) design matrix
    A_sp    = sparse(numel(nFx)*n_img,length(Z));
    is_meas = false(size(A_sp,1),1);
    for i = 1 : n_img % Loop through angles of illumination
        pupDom(:,:,i) = pupDom(:,:,i).*(nFx.^2 + nFy.^2 >= 0.05); % 0 out low-frequencies

        MsCell{i}   = find(pupDom(:,:,i));

        [polyTrans,AP,AM] = getPolyTransPeriodogram( maxDeg, U(i), V(i), Z );
        % Design matrix
        A = BMat(MsCell{i},:)*polyTrans;
        nonzeroDesign{i} = (sum(abs(A)) ~= 0);

        A_sp(MsCell{i}+numel(nFx)*(i-1),nonzeroDesign{i}) = A(:,nonzeroDesign{i});
        is_meas(MsCell{i}+numel(nFx)*(i-1)) = true;
    end
end