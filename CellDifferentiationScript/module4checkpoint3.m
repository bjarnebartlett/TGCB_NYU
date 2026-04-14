%% ============================================================
% HEMATOPOIETIC DIFFERENTIATION (MATCHES IMAGE EXACTLY)
%% ============================================================

modelName = 'Hematopoiesis_Model';

%% 1. Cleanup
if bdIsLoaded(modelName)
    close_system(modelName,0);
end

new_system(modelName);
open_system(modelName);
rt = sfroot;

%% 2. Add Chart
add_block('sflib/Chart',[modelName '/Hematopoiesis_Logic']);
chart = rt.find('-isa','Stateflow.Chart','Path',...
    [modelName '/Hematopoiesis_Logic']);

%% ============================================================
%% 3. STATE POSITIONS (Manually spaced to match diagram)
%% ============================================================

% Core lineage
HSC = makeState(chart,'Hematopoietic_Stem_Cell',[100 300 230 60]);
MSC = makeState(chart,'Multipotent_Stem_Cell',[450 300 240 60]);

LP  = makeState(chart,'Lymphoid_Progenitor',[850 150 220 60]);
MP  = makeState(chart,'Myeloid_Progenitor',[850 450 220 60]);

% Lymphoid outputs
NK = makeState(chart,'NK_Cell',[1200 50 170 50]);
T  = makeState(chart,'T_Cell',[1200 130 170 50]);
B  = makeState(chart,'B_Cell',[1200 210 170 50]);

% Myeloid outputs
Neut = makeState(chart,'Neutrophil',[1200 350 180 50]);
Baso = makeState(chart,'Basophil',[1200 420 180 50]);
Eos  = makeState(chart,'Eosinophil',[1200 490 180 50]);
Mono = makeState(chart,'Monocyte',[1200 560 180 50]);
RBC  = makeState(chart,'Red_Blood_Cell',[1200 630 200 50]);
Plat = makeState(chart,'Platelet',[1200 700 180 50]);

%% ============================================================
%% 4. TRANSITIONS (Exact biological flow)
%% ============================================================

% Default → HSC
dt = Stateflow.Transition(chart);
dt.Destination = HSC;

% HSC → Multipotent Stem Cell
addTransition(HSC,MSC,'[stem_activation == 1]',3,9);

% Multipotent → Lymphoid Progenitor
addTransition(MSC,LP,'[lymphoid_signal == 1]',2,9);

% Multipotent → Myeloid Progenitor
addTransition(MSC,MP,'[myeloid_signal == 1]',4,9);

% Lymphoid differentiation
addTransition(LP,NK,'[IL15 > 0.5]',3,9);
addTransition(LP,T,'[Notch_signal == 1]',3,9);
addTransition(LP,B,'[IL7 > 0.5]',3,9);

% Myeloid differentiation
addTransition(MP,Neut,'[GCSF > 0.5]',3,9);
addTransition(MP,Baso,'[IL3 > 0.5]',3,9);
addTransition(MP,Eos,'[IL5 > 0.5]',3,9);
addTransition(MP,Mono,'[MCSF > 0.5]',3,9);
addTransition(MP,RBC,'[EPO > 0.5]',3,9);
addTransition(MP,Plat,'[TPO > 0.5]',3,9);

%% ============================================================
%% 5. INPUT SIGNALS
%% ============================================================

inputs = {'stem_activation','lymphoid_signal','myeloid_signal',...
          'IL15','IL7','Notch_signal','GCSF','IL3','IL5',...
          'MCSF','EPO','TPO'};

for i = 1:length(inputs)
    d = Stateflow.Data(chart);
    d.Name = inputs{i};
    d.Scope = 'Input';
end

%% ============================================================
%% 6. Fit View
%% ============================================================

chart.Editor.open;
pause(0.5);
sfroot.getCurrentEditor.fitToView;

save_system(modelName);
disp('Biologically Correct Hematopoiesis Model Created.');

%% ============================================================
%% Helper Functions
%% ============================================================

function s = makeState(parent,name,pos)
    s = Stateflow.State(parent);
    s.Name = name;
    s.Position = pos;
end

function addTransition(src,dst,label,so,do)
    t = Stateflow.Transition(src.Chart);
    t.Source = src;
    t.Destination = dst;
    t.LabelString = label;
    t.SourceOclock = so;
    t.DestinationOclock = do;
end