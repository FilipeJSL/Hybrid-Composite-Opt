function [Fn] = MapVariables(Fn, AdmissibleSet)

%Maps the GA design variables into the actual set of fibres available in
%each of the fibre groups chosen in OptModule.m
Fn(1)=AdmissibleSet{1,1}(Fn(1));
Fn(2)=AdmissibleSet{1,2}(Fn(2));

end

