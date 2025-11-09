%%
clear cutoff delta e_field endtime fm_times Fs fs1 fs2 keithley loc pk time

%fs1 = 'poly_drop_K.txt';
%fs2 = 'poly_drop_S.dat';
%fs1 = 'uncharged_glass_drop_K.txt';
%fs2 = 'uncharged_glass_drop_S.dat';
%fs1 = 'mgs1_drop_KM.txt';
%fs2 = 'mgs1_drop_S.dat';
fs1 = 'mgs1c_drop_KM.txt';
fs2 = 'mgs1c_drop_S.dat';


[keithley,time,delta,Fs,cutoff,endtime,M] = keithley_import(fs1,20);
%%
[e_field,fm_times,pk,loc] = fm_import(fs2,cutoff,endtime);

%%
if fs1 == "poly_drop_K.txt"
    loc = 2.7;
elseif fs1 == "uncharged_glass_drop_K.txt"
    loc = 4.6;
elseif fs1 == "mgs1_drop_KM.txt"
    loc = 3.7;
elseif fs1 == "mgs1c_drop_KM.txt"
    loc = 7.8;
else
    loc = 0;
end


time = time - loc;
fm_times = fm_times - loc;
%%
figure
yyaxis left
plot(time,keithley) 
ylabel('Charge [pC]')
yyaxis right
plot(fm_times,e_field,':')
hold on
lowess = smooth(fm_times,e_field,30,'lowess');
plot(fm_times,lowess,'-','LineWidth',3);
ylabel('Electric Field [V/m]')
xlabel('Time [s]')
%%
figure
plot(time,delta)
%%
function[down_keithley,down_time,delta_keithley,Fs,cutoff,endtime,M] = keithley_import(fs1,range)
    opts = delimitedTextImportOptions("NumVariables", 7);
    opts.DataLines = [1, Inf];
    opts.Delimiter = ",";
    opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string"];
    temp = readmatrix(fs1,opts);
    range_fac = range/2;
    keithley = str2double(temp(4:end,7)) * range_fac;
    keithley = keithley * (28e-12/1e-12); %capacitance is 28e-12 and convert to pC
    Fs = str2double(temp(2,2));
    M = str2double(temp(4:end,1:6));
    time = transpose(linspace(0,length(keithley)/Fs, length(keithley)));
    cutoff = datetime(temp(1,1),'format','dd/MM/uuuu HH:mm:ss.SSSSSS');
    endtime = time(end);
    max_keithley = 10;
    block = Fs / max_keithley;
    new_l = length(keithley)/block;
    down_keithley = NaN(new_l,1);
    down_time = linspace(0,length(keithley)/Fs,new_l);
    delta_keithley = NaN(new_l,1);
    for i = 1:new_l
        down_keithley(i) = mean(keithley(1+block*(i-1):block*i));
    end
    for i = 2:new_l
        delta_keithley(i) = down_keithley(i) - down_keithley(i-1);
    end
end

function [e_field,fm_times,pk,loc] = fm_import(fs2,cutoff,endtime)
    opts = delimitedTextImportOptions("NumVariables", 8);
    opts.DataLines = [5, Inf];
    opts.Delimiter = ",";
    opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string"];
    temp = readmatrix(fs2,opts);
    e_field = str2double(temp(:,3));
    times = temp(:,1);
    for i = 1:length(times)
        if times{i}(end-1) ~= '.'
            times{i} = strcat(times{i},'.0');
        end
    end
    fm_times = datetime(times,'format','dd/MM/uuuu HH:mm:ss.S');
    e_field = e_field(fm_times>cutoff);
    fm_times = fm_times(fm_times>cutoff);

    fm_times = seconds(datetime(fm_times) - datetime(fm_times(1)));

    e_field = e_field(fm_times<endtime);
    fm_times = fm_times(fm_times<endtime);

    [pk,loc] = findpeaks(e_field,'MinPeakHeight',0.5);
    if isempty(pk)
        [pk,loc] = findpeaks(-e_field,'MinPeakHeight',0.5);
    end
    if length(pk) > 1
        [pk,loc] = findpeaks(e_field,'MinPeakHeight',1);
        if isempty(pk)
            [pk,loc] = findpeaks(-e_field,'MinPeakHeight',1);
        end
    end
    loc = fm_times(loc);
end