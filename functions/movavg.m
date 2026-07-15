function s2 = movavg(s1, Sr)

% s1 has to be along  by column signals
% Sr is the Sampling Rate

[N, M] = size(s1);

% initializing the output
s2 = zeros(N, M);

% updating the output's values
for i = 1:N
    in = max(1, i-Sr);
    s2(i,:) = mean(s1(in:i,:));
end
