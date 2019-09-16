% Wrapping function, which takes a scalar, vector, or matrix, and
% reassigns all of its elements to fall into a given closed interval, e.g.:
%
% interval_wrap((-10):10, [-5 5]) =
% [0 1 2 3 4 -5 -4 -3 -2 -1 0 1 2 3 4 5 -4 -3 -2 -1 0]
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% IN:   Data         =  data (scalar, vector, or matrix)
% IN:   Interval     =  two-element vector of the interval in which Data is to
%                       be wrapped; if a data element falls within (or on a
%                       boundary of) the interval, it is unchanged; otherwise,
%                       it is brought up or down in steps of
%                       Interval(2) - Interval(1).
%                       Note that Interval(2) > Interval(1) is required.
%
% OUT:  Data_wrapped =  The result of wrapping Data as described above.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%

function Data_wrapped = interval_wrap(Data, Interval)

if Interval(2) <= Interval(1)
    error('Interval(2) > Interval(1) is required')
else
end

if numel(size(Data)) > 2
    error('I have not coded this to act on 3+Dimensional matrices.')
else
end

D = Interval(2) - Interval(1);

% % This works, at least for vectors, but it is very slow for matrices:
% while sum(Data<=Interval(2) & Data>=Interval(1)) < numel(Data)
%     Data(Data<Interval(1)) = Data(Data<Interval(1)) + D;
%     Data(Data>Interval(2)) = Data(Data>Interval(2)) - D;
% end

if isvector(Data) % true for scalars, rows, and columns
    for i=1:length(Data)
        if Data(i) < Interval(1) || Data(i) > Interval(2)
            Data(i) = Data(i) + D*round((mean(Interval) - Data(i))/D);
        else
        end
    end
else
    for i=1:size(Data,1)
        for j=1:size(Data,2)
            if Data(i,j) < Interval(1) || Data(i,j) > Interval(2)
                Data(i,j) = Data(i,j) + D*round((mean(Interval) - Data(i,j))/D);
            else
            end
        end
    end
end

Data_wrapped = Data;

end
