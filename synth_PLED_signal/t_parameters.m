CUR_DIR=mfilename('fullpath');
dir_st=(dir([CUR_DIR '.m']));
CUR_DIR=CUR_DIR(1:end-(length(dir_st.name)-2));

DATA_DIR=[CUR_DIR '..' filesep 'data' filesep];
RES_DIR=[DATA_DIR 'results' filesep];
FNAME_STUB='synth_PLEDs_';
IMP_STR='I';

DUFFING_WAVEFORMS_FILE=[DATA_DIR 'duffing_waveforms.mat'];

FONT_TYPE='Arial';
FONT_SIZE=13;




