function ml_report_writer(outputfile, pvalue, ts, both, ... 
    featmat1, featmat2, setname1, setname2, slfnames, slfids, ... 
    slfdescrp, conf, id, succ,localhost,mvtestmeth)

%Obsolete

error(tz_genmsg('of','ml_report_writer','ml_simec_report'));

nfeats1=size(featmat1,2);
nfeats2=size(featmat2,2);

nsample1=size(featmat1,1);
nsample2=size(featmat2,1);

fid = fopen(outputfile, 'w');
fprintf(fid, '<HTML><HEAD><TITLE>SImEC Results Report</TITLE></HEAD>\n');
fprintf(fid, '<BODY aLink=#0077ff bgColor=#ffffff link=#118811 text=#000000 vLink=#ff0000>\n');  

fprintf(fid, '<CENTER><B><U>SImEC Results Report'); 
fprintf(fid, '.</U></B></CENTER>\n<BR>\n\n'); 

if (id ~= -1) 
    st_time = datestr(id/10000); 
    fprintf(fid, [' for SImEC session started on&nbsp;' st_time ' GMT']); 
end

if succ<=0
    fprintf(fid,'<P>The two sets <U>%s</U> and <U>%s</U> have been selected for comparison.</P>',...
        setname1,setname2);
    
    fprintf(fid,'<H3>ERROR:</H3>');
    
    if (nsample1+nsample2)<nfeats1
        succ=-3;
    end
    if nfeats1 ~= nfeats2
        succ=-2;
    end
    
    switch succ
    case 0
        fprintf(fid,'<P>Comparison can not be not completed because of unknown errors.</P>');
    case -1
        fprintf(fid,'<P>Some errors happend in individual feature comparison.</P>');
    case -2
        fprintf(fid,'<P>The total number of images is less than the number of features.</P>');
    case -3
         fprintf(fid,'<P>The numbers of features in the two sets do not equal.</P>');   
    end
else
    fprintf(fid,'<H3>Summary</H3>');
    fprintf(fid,['<P><B>Multivariate testing method: </B>Two-sample ', mvtestmeth, '</P>']);
    if pvalue<=conf
        fprintf(fid,'<P>The two sets <U>%s</U> and <U>%s</U> are different at the conficence level <B>%s</B>. ',setname1,setname2,num2str(1-conf));
    else
        fprintf(fid,'<P>The two sets <U>%s</U> and <U>%s</U> are NOT concluded to be different at the confidence level <B>%s</B>. ',setname1,setname2,num2str(1-conf));
    end
    
    fprintf(fid,'<BR>The p-value is <B>%3.3f</B>.',pvalue);
    
    if pvalue<=conf
        fprintf(fid,'<H3>Features ranked by p-values:</H3>');
        fprintf(fid,'<P><B>Univariate testing method: </B>Two-sample t-test</P>');
        fprintf(fid, '<TABLE align=center cellpadding=2><TBODY>\n');
        fprintf(fid, '<TR>\n');
        fprintf(fid,'%s%s%s%s\n','<TD><B><U>Feature ID</B></U></TD>','<TD><B><U>Feature Name</B></U></TD>','<TD><B><U>Feature Description</B></U></TD>', '<TD><B><U>p-value</B></U></TD>');
        fprintf(fid, '</TR>\n');
        
        both(:,both(2,:)<0.5)=[];
        
        for i = 1:size(both,2)
            fprintf(fid, '<TR>\n');
            slfid=slfids{both(1,i)};
            if str2num(slfid)>23.5 & str2num(slfid)<72.5
                idlink='Z';
            else
                idlink=slfid;
            end
            fprintf(fid,'%s%12.10f%s\n',['<TD>',slfid, '</TD><TD>','<A href="',localhost,'SLF/features.html#', idlink, '">',slfnames{both(1,i)},'</A></TD><TD>',slfdescrp{both(1,i)},'</TD><TD>'],both(3,i), '</TD>');
            %end
            fprintf(fid, '</TR>\n');
        end
        
        fprintf(fid, '</TBODY></TABLE>');
    end
end
fprintf(fid, '</BODY></HTML>\n');  
status = fclose(fid); 
