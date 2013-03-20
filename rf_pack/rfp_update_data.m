function out = rfp_update_data(S)

for isession=1:length(S.session)
    session = S.session(isession);
    norients = length(session.orientations);
    for icell=1:length(session.cells)
        c = session.cells(icell);
        joint_trains = {};
        for iorient=1:norients
            s_trains = c.spike_trains(iorient, :);
            tmp_train = [];
            for itrain=1:length(s_trains)
                train = s_trains{itrain};
                tmp_train = [tmp_train; train]; %#ok
            end
            joint_trains{end+1} = sort(tmp_train); %#ok            
        end
        S.session(isession).cells(icell).joint_trains = joint_trains;
    end
end

out = S;