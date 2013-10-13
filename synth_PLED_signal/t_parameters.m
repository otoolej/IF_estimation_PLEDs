DATA_DIR='~/deusto/software/EEG/PLEDs/data/synthetic_PLEDs/';
RES_DIR=[DATA_DIR 'results/'];
FNAME_STUB='synth_PLEDs_';
IMP_STR='I';


CUR_DIR=mfilename('fullpath');
dir_st=(dir([CUR_DIR '.m']));
CUR_DIR=CUR_DIR(1:end-(length(dir_st.name)-2));

DATA_DIR=[CUR_DIR '..' filesep 'data' filesep];

DUFFING_WAVEFORMS_FILE=[DATA_DIR 'duffing_waveforms.mat'];


FONT_TYPE='Times-Roman';
FONT_SIZE=14;

