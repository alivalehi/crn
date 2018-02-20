function update_analysis(fname, save_active)
%cd data
readf = dir(fname);
nSIM = length (readf);

for i = 1: nSIM
    filename = readf(i).name;
    load(filename);
    fprintf('File:%s  N:%d H:%d PER:%1.3f \n', filename, sim.N, sim.H, sim.Exp_PER);  
    
    sim = fanalyzesim(sim);
    if (save_active)
        fprintf('File analysis updated and saved.\n\n');
        save(filename, 'sim');
    else
        fprintf('File analysis performed.\n\n');
    end
    
end
%cd ..
disp('end');
pause;
end