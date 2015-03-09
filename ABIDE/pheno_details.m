close all
clear all
fID=fopen('Phenotypic_V1_0b.csv');
% SITE_ID,SUB_ID,DX_GROUP,DSM_IV_TR,AGE_AT_SCAN,SEX,HANDEDNESS_CATEGORY,HANDEDNESS_SCORES,FIQ,VIQ,PIQ,FIQ_TEST_TYPE,VIQ_TEST_TYPE,PIQ_TEST_TYPE,ADI_R_SOCIAL_TOTAL_A,ADI_R_VERBAL_TOTAL_BV,ADI_RRB_TOTAL_C,ADI_R_ONSET_TOTAL_D,ADI_R_RSRCH_RELIABLE,ADOS_MODULE,ADOS_TOTAL,ADOS_COMM,ADOS_SOCIAL,ADOS_STEREO_BEHAV,ADOS_RSRCH_RELIABLE,ADOS_GOTHAM_SOCAFFECT,ADOS_GOTHAM_RRB,ADOS_GOTHAM_TOTAL,ADOS_GOTHAM_SEVERITY,SRS_VERSION,SRS_RAW_TOTAL,SRS_AWARENESS,SRS_COGNITION,SRS_COMMUNICATION,SRS_MOTIVATION,SRS_MANNERISMS,SCQ_TOTAL,AQ_TOTAL,COMORBIDITY,CURRENT_MED_STATUS,MEDICATION_NAME,OFF_STIMULANTS_AT_SCAN,VINELAND_RECEPTIVE_V_SCALED,VINELAND_EXPRESSIVE_V_SCALED,VINELAND_WRITTEN_V_SCALED,VINELAND_COMMUNICATION_STANDARD,VINELAND_PERSONAL_V_SCALED,VINELAND_DOMESTIC_V_SCALED,VINELAND_COMMUNITY_V_SCALED,VINELAND_DAILYLVNG_STANDARD,VINELAND_INTERPERSONAL_V_SCALED,VINELAND_PLAY_V_SCALED,VINELAND_COPING_V_SCALED,VINELAND_SOCIAL_STANDARD,VINELAND_SUM_SCORES,VINELAND_ABC_STANDARD,VINELAND_INFORMANT,WISC_IV_VCI,WISC_IV_PRI,WISC_IV_WMI,WISC_IV_PSI,WISC_IV_SIM_SCALED,WISC_IV_VOCAB_SCALED,WISC_IV_INFO_SCALED,WISC_IV_BLK_DSN_SCALED,WISC_IV_PIC_CON_SCALED,WISC_IV_MATRIX_SCALED,WISC_IV_DIGIT_SPAN_SCALED,WISC_IV_LET_NUM_SCALED,WISC_IV_CODING_SCALED,WISC_IV_SYM_SCALED,EYE_STATUS_AT_SCAN,AGE_AT_MPRAGE,BMI
% we keep columns 2,3,5,6,9,15,16,17,18,19,21,22,23,24,25
% that is:
% SUB_ID,DX_GROUP,AGE_AT_SCAN,SEX,FIQ,ADI_R_SOCIAL_TOTAL_A,ADI_R_VERBAL_TOTAL_BV,ADI_RRB_TOTAL_C,ADI_R_ONSET_TOTAL_D,ADI_R_RSRCH_RELIABLE,ADOS_TOTAL,ADOS_COMM,ADOS_SOCIAL,ADOS_STEREO_BEHAV,ADOS_RSRCH_RELIABLE

keepcols=[2,3,5,6,9,15,16,17,18,19,21,22,23,24,25];
str='';
for k=1:max(keepcols);
    if(ismember(k,keepcols))
        str=[str ' %s'];
    else
        str=[str ' %*s']; % matlab way of skipping a column
    end
end
a=textscan(fID,[str ' %*[^\n]'],'Delimiter',',','HeaderLines',1); % now a has 15 cells for each column
fclose(fID);

% loading the sites
fID=fopen('Phenotypic_V1_0b.csv');
keepcols=[1];
str='';
for k=1:max(keepcols);
    if(ismember(k,keepcols))
        str=[str ' %s'];
    else
        str=[str ' %*s']; % matlab way of skipping a column
    end
end

sites=textscan(fID,[str ' %*[^\n]'],'Delimiter',',','HeaderLines',1); % now a has 15 cells for each column


% convert them into a matrix

phenodata=zeros(length(a{1}),length(a));
for c=1:length(a);
    for r=1:size(phenodata,1)
        temp=str2num(a{c}{r});
        if(isempty(temp))
            temp=-9999;
        end
        phenodata(r,c)=temp;
    end
end

% filter the matrix
ASD=phenodata(:,2)==1; % keep patients
sex=phenodata(:,4)==1; % keep males
iq=phenodata(:,5)>=100;	% keep high functioning
score_columns=6:15; % 6:10 is ADI-R 11:15 is ADOS
adi_ados=phenodata(:,score_columns)>0; % keep subjects with ADI-R and ADOS

goodsites={'CALTECH','CMU','PITT','UM_1','UM_2'}; % same sites as in http://www.nature.com/neuro/journal/v18/n2/full/nn.3919.html
for s=1:length(goodsites)
    keepsites(:,s)=strcmp(sites{:},goodsites{s});
end

% for compatibility with old matlab, we need to convert logical variables to single (or double)
list=find(single(ASD).*single(sex).*single(iq).*prod(single(adi_ados),2).*sum(single(keepsites),2)); % the final rows corresponding to subjects that have good data
scores=phenodata(list,score_columns);
score_labels={'ADI_R_SOCIAL_TOTAL_A', 'ADI_R_VERBAL_TOTAL_BV', 'ADI_RRB_TOTAL_C', 'ADI_R_ONSET_TOTAL_D', 'ADI_R_RSRCH_RELIABLE', 'ADOS_TOTAL', 'ADOS_COMM', 'ADOS_SOCIAL', 'ADOS_STEREO_BEHAV', 'ADOS_RSRCH_RELIABLE'};
%save scores scores score_labels;
