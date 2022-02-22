function envout=envsm(trin,t,op_length,trip)

%  envsm=smooth(trin,t,op_length)
%
% envsm computes a smooth amplitude envelope using the following method:
% 
%   1) Compute Hilbert envelope of trin
%   2) Convolve envelope with triangular smoother of half-length
%       op_length
%
% The program is based on aec.m written by Gary Margrave
%
% trin= input trace
% t= time coordinate vector for trin
% op_length= half-length of triangular smoother in same units as t
%
% enhvsm = smoothed output envelope
%
% D. Eaton, October 2010
   
% set defaults
% double the operator length
 op2=op_length*2;
% form new trace padded to a power of 2
 trinew=padpow2(trin,0);
% compute the envelope
 env=abs(hilbm(trinew));
 env=env(1:length(trin));
% guard against accidental transpose
           if size(env)~=size(trin), env=env.';end
% compute the smoothed envelope
nop=round(op2/(t(2)-t(1)))+1;
iimax = floor(nop/2);
for ii = 1:iimax
   triang(ii) = ii*2/nop;
   triang(nop-ii+1) = triang(ii);
end
triang(iimax+1) = 1.0;

triang = triang/length(triang);

 envout=conv(env,triang);
% grab the central length(trin) samples
 envout=envout(round(nop/2):length(trin)+round(nop/2)-1);









