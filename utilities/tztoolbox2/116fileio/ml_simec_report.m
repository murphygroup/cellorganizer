function ml_simec_report(outputfile, pvalue, both, ... 
    featmat1, featmat2, setname1, setname2, ... 
    featuresetname, conf, succ,localhost,mvtestmeth,featids)
%ML_SIMEC_REPORT generats an html page for simec results
%   ML_SIMEC_REPORT(OUTPUTFILE,PVALUE,BOTH,FEATMAT1,FEATMAT2,
%       SETNAME1,SETNAME2,FEATURESETNAME,
%	CONF,SUCC,LOCALHOST,MVTESTMETH)
%   creates a file with the full path is OUTPUTFILE. This file reports
%   results of comparison. PVALUE is the p-value of comparison.
%   BOTH is a result matrix of univariate testings and has two rows. 
%   The first row contains pvalues and the second row contains feature 
%   indieces. FEATMAT1 and FEATMAT2 are feature matrices. SETNAME1 and 
%   SETNAME2 are set names. LOCALHOST is the server for feature description.
%   MVTESTMETH is the multivariate testing method.

if ~exist('featids','var')
    featids = [];
end

nfeats1=size(featmat1,2);
nfeats2=size(featmat2,2);

nsample1=size(featmat1,1);
nsample2=size(featmat2,1);

fid = fopen(outputfile, 'w');
fprintf(fid, '<HTML><HEAD><TITLE>SImEC Results Report</TITLE></HEAD>\n');
fprintf(fid, ['<BODY aLink=#0077ff bgColor=#ffffff link=#118811' ...
    'text=#000000 vLink=#ff0000>\n']);  

fprintf(fid, '<CENTER><B><U>SImEC Results Report'); 
fprintf(fid, '.</U></B></CENTER>\n<BR>\n\n'); 

if succ<=0
    fprintf(fid,['<P>The two sets <U>%s</U> and <U>%s</U> ' ...
        'have been selected for comparison.</P>'],setname1,setname2);
    
    fprintf(fid,'<H3>ERROR:</H3>');
    
    if (nsample1+nsample2)<nfeats1
        succ=-3;
    end
    if nfeats1 ~= nfeats2
        succ=-2;
    end
    
    switch succ
    case 0
        fprintf(fid,['<P>Comparison can not be not completed' ...
            'because of unknown errors.</P>']);
    case -1
        fprintf(fid,'<P>Some errors happend in individual feature comparison.</P>');
    case -2
        fprintf(fid,'<P>The total number of images is less than the number of features.</P>');
    case -3
         fprintf(fid,'<P>The numbers of features in the two sets do not equal.</P>');   
    end
else
    fprintf(fid,'<H3>Summary</H3>');
    fprintf(fid,['<P><B>Multivariate testing method: </B>Two-sample ', ...
        mvtestmeth, '</P>']);
    fprintf(fid,['<P><B>Feature Set: <B>']);
    if isempty(localhost)
        fprintf(fid,'</P>');
    else
        fprintf(fid,['<A href="',localhost,'SLF/" target="_blank">', ...
            featuresetname, '</A></P>']);
    end
    
    if pvalue<=conf
        fprintf(fid,['<P>The two sets <U>%s</U> and <U>%s</U> ' ...
            'are different at the confidence level <B>%s</B>. '], ...
            setname1,setname2,num2str(1-conf));
    else
        fprintf(fid,['<P>The two sets <U>%s</U> and <U>%s</U> are NOT ' ...
            'concluded to be different at the confidence level <B>%s</B>. '], ...
            setname1,setname2,num2str(1-conf));
    end
    
    fprintf(fid,['<BR>The p-value is <B>' num2str(pvalue) '</B>.']);
    
    if pvalue<=conf
        fprintf(fid,'<H3>Features ranked by p-values:</H3>');
        fprintf(fid,'<P><B>Univariate testing method: </B>Two-sample t-test</P>');

        both(:,both(2,:)<0.5)=[];
        [featdesp,featsetnames,allfeatids] = ml_featinfo;
        selfeatset = strmatch(featuresetname,featsetnames,'exact');
        
        if ~isempty(featids)
            [tmp,featididx] = ismember(sort(featids), ...
                                       allfeatids{selfeatset});
        else
            [tmp,featididx] = sort(allfeatids{selfeatset});
        end
        
        s3 = featididx(both(1,:));
        fprintf(fid,ml_showfeatinfo_html( ...
            featuresetname,s3,'p-values',both(3,:)));
    end
end
fprintf(fid, '</BODY></HTML>\n');  
status = fclose(fid); 
