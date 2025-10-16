function rl2 = rl2Err(IH,data)
%rL2Err: Calculates relative l2 error percentage between simulated and
%           observed bacterial count.
%   Inputs:
%       OD      <- simulated bacterial count
%       data    <- observed bacterial count
%   Outputs:
%       rl2Err  <- relative l2 error percentage

diff = IH - data';
rl2 = sqrt(sum(sum(diff.^2))/sum(sum(data.^2))) * 100;

end