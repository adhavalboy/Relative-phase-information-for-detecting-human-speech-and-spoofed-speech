function [feature_cmv]= Relative_phase_shift(x,fs,freqScale)
x = filter([1 -0.97],1,x);
[srh_f0,srh_vuv,srh_vuvc,srh_time] = pitch_srh(x,fs,70,700,10);

fos = [srh_time',srh_f0'];

opt = sin_analysis();
opt.fharmonic  = true; % Use harmonic model
opt.use_ls = false; % Use Peak Picking
frames = sin_analysis(x, fs, fos,opt);
opt.pd_vtf_rm    = false; % Remove the VTF phase from the inst. phase

opt.dc_phase     = 0;    % if empty, does not change it; otherwise, set the given value
opt.polarity_inv = true;% For vizualisation purpose, the phase of the
                         % inverted signal might be more convenient.
                         % (applied after DC's phase is set using dc_phase)
opt.pd_method    = 2;    % 1:Phase Distortion (PD) [1-3]
                         % 2:Relative Phase Shift (RPS) [4-5]
opt.harm2freq    = true;% Convert the harmonic values on a hertz scale

opt.dftlen    = 512;    % The DFT length used for envelope estimation

opt.usemex    = false; % Use mex fn, faster but use linear interpolation
opt.debug = 0;



[RPS, AE, opt] = phase_rpspd(frames, fs, opt);
N=4096; t=0:1/fs:(length(x)-1)/fs;
f=0:fs/N:fs/2; Fmin = 300; Fmax = fs/2;
if strcmp(freqScale,'mel')
    low = 1125*log(1+(Fmin/700));
    high = 1125*log(1+(Fmax/700));
    frqsc = linspace(low,high,70);
    hz = 700*(exp(frqsc/1125)-1);
elseif strcmp(freqScale,'linear')
    low = Fmin;
    high = Fmax;
    hz = linspace(low,high,70);
    frqsc = hz;
elseif strcmp(freqScale,'log10')
    low = log10(Fmin);
    high = log10(Fmax);
    frqsc = linspace(low,high,70);
    hz = 10.^frqsc;
elseif strcmp(freqScale,'erb')
    low = frq2erb(Fmin);
    high = frq2erb(Fmax);
    frqsc = linspace(low,high,70);
    hz = erb2frq(frqsc);
elseif strcmp(freqScale,'bark')
    low = frq2bark(Fmin);
    high = frq2bark(Fmax);
    frqsc = linspace(low,high,70);
    hz = bark2frq(frqsc);
end
fil_res = floor((N+1)*hz/fs);
for k=1:length(f)
        for m=1:length(frqsc)-2
            if(k<fil_res(m))
                h(m,k)=0;
            elseif(fil_res(m)<=k && k<=fil_res(m+1))
                h(m,k)=(k-fil_res(m))/(fil_res(m+1)-fil_res(m));
            elseif(fil_res(m+1)<=k && k<=fil_res(m+2))
                h(m,k)=(fil_res(m+2)-k)/(fil_res(m+2)-fil_res(m+1));
            elseif(k>fil_res(m+2))
                h(m,k)=0;
            end
        end
end
for k=1:size(RPS,1)
    uw(:,k) = unwrap(RPS(k,:));
    RPS1(k,:) = uw(:,k)'*h';
    DPE(k) = mean(diff([0 RPS1(k,:)]));
end
RP1 = RPS1';
tem = [DPE;RP1];
temp = dct2(tem);
static = temp(1:60,:);
del = deltas(static);
deldel = deltas(del);
feature = [static;del;deldel];
feature_cmv = cmvn(feature,'false');