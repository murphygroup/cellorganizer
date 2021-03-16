%Thhis script is to test how cvs is working
a=1;
b=2;
c=3;
d=4;

tz_save([tz_getdir('data') filesep 'cvstest.mat'],{a,b,c}, ...
    {'a','b','c'},'tz_cvstest','for testing cvs');