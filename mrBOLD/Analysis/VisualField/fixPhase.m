function view = fixPhase(view, scanNum, sliceNum, shift, scale)%% function view = fixPhase(view, [scanNum], [sliceNum], [shift], [scale])%%% Adds specified offset and scale to the phases of the specified coranal.% Will prompt for offset & scale if not given of empty.% The current slice will be used if sliceNum not given or empty.%% HISTORY:%   2002.03.08 RFD (bob@white.stanford.edu) wrote it.global dataTYPES;if(~exist('scanNum','var') | isempty(scanNum))    scanNum = getCurScan(view);endif(~exist('sliceNum','var') | isempty(sliceNum))    sliceNum = viewGet(view, 'Current Slice');endoldShift = dataTYPES(view.curDataType).atlasParams(scanNum).phaseShift(sliceNum);oldScale = dataTYPES(view.curDataType).atlasParams(scanNum).phaseScale(sliceNum);oldPh = view.ph;% undo the old shift & scale so that we can start anewview.ph{scanNum}(:,:,sliceNum) = mod((view.ph{scanNum}(:,:,sliceNum) - oldShift)/oldScale, 2*pi);%refreshView(view);if(~exist('shift','var')) shift = []; endif(~exist('scale','var')) scale = []; endif(isempty(shift) | isempty(scale))    if(isempty(shift))        def{1} = num2str(oldShift);    else        def{1} = num2str(shift);    end    if(isempty(scale))        def{2} = num2str(oldScale);    else        def{2} = num2str(scale);    end        ok = 0;    while(~ok)        answer = inputdlg({'Phase shift:','Scale Factor:'}, 'Fix Phase', 1, def);        if(isempty(answer))            % Cancel all changes            view.ph = oldPh;            refreshView(view);            ok = 1;            return;        end        shift = eval(answer{1});        scale = eval(answer{2});        view.ph{scanNum}(:,:,sliceNum) = view.ph{scanNum}(:,:,sliceNum) * scale + shift;        view.ph{scanNum}(:,:,sliceNum) = phaseStandardize(view.ph{scanNum}(:,:,sliceNum));        view.ph{scanNum}(:,:,sliceNum) = mod(view.ph{scanNum}(:,:,sliceNum), 2*pi);        refreshView(view);        ans = questdlg('Is this OK?', ...            'Confirm Phase Adjustment', ...            'Yes','No- try again','Cancel all changes','Yes');                switch ans,        case 'Yes',             saveCorAnal(view);            ok = 1;        case 'No- try again',            refreshView(view);            def = answer;            ok = 0;        case 'Cancel all changes',            view.ph = oldPh;            refreshView(view);            ok = -1;            return;        end % switch    endelse    view.ph{scanNum}(:,:,sliceNum) = view.ph{scanNum}(:,:,sliceNum) * scale + shift;    view.ph{scanNum}(:,:,sliceNum) = phaseStandardize(view.ph{scanNum}(:,:,sliceNum));    view.ph{scanNum}(:,:,sliceNum) = mod(view.ph{scanNum}(:,:,sliceNum), 2*pi);    refreshView(view);    saveCorAnal(view);enddataTYPES(view.curDataType).atlasParams(scanNum).phaseShift(sliceNum) = shift;dataTYPES(view.curDataType).atlasParams(scanNum).phaseScale(sliceNum) = scale;saveSession;return;function ph = phaseStandardize(ph)% wrap neagative phases back up to the 2pi endl = ph<0;if(~isempty(l))    ph(l) = ph(l)+2*pi;endreturn;