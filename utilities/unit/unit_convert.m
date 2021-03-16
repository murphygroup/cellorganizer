function varargout = unit_convert(varargin)
%UNIT_CONVERT converts between units, also checks dimensions
%   UNIT_CONVERT(str1,str2) displays the conversion factor between str1 and str2.
%   examples: UNIT_CONVERT('N','kg.m/s2') or UNIT_CONVERT('MW','J/s')
%   note the units can be compound terms with powers and prefixes
%
%   UNIT_CONVERT(str1,str2,num) converts the specified number of units
%
%   f        = UNIT_CONVERT(...) returns the conversion as an output
%   [f,resu] = UNIT_CONVERT(...) returns additional information about the units
%
%   UNIT_CONVERT('report') displays a report of the code:
%                  units      (e.g. meter)
%                  unit kinds (e.g. length)
%                  prefixes   (e.g. M for mega=10^6)
%                  options    (e.g. 'abstemp' for absolute temperatures)
%
%   other options: list, report, show, verb, detail, eng, abstemp, showfrac
%
%EXAMPLES
%   unit('kg','lb')       %converts from kilos to pounds
%   unit('atm','lbf/in2') %converts from standard atmospheres to psi.
%   unit('lbf/lb','m/s2') %acceleration due to gravity in m/s^2
%   unit('lbf-in','N.m')  %moment conversion
%   unit('cc','m3',3.2)   %convert 3.2 cc to m3
%   unit('N.m','kg*m2/s^2')
%   unit('W','kW.hr/day') %kiloWatt-hours per day
%   unit('rpm','Hz')      %Hz is rev/second, not just 1/second
%
%UNITS
%   unit('report') displays all units recognised by this function
%   unit('list')   displays just their names
%
%   these may be compound units, e.g. N.m or lb/s
%   terms must be separated by one of . * / or space
%
%   each term may have a power, e.g m^3 or m3. the ^ is optional.
%
%   each term may have a scaling prefix, e.g. kW for kiloWatt
%   powers are applied to prefixes: e.g. km3 is km x km x km
%
%   divisors are applied to just the following term
%   kg/m*s2 is kg^1 m^-1 s^2 : note how the power of "s" is positive
%
%CHECKS
%   UNIT_CONVERT checks that input/output units have identical dimensions
%   e.g. that both are velocities, or mass/time, etc.
%   this is done by reducing each unit to basic "kinds" of unit, e.g length
%   unit('report') will display all the kinds of unit
%
%CUSTOM UNITS
%   it is possible to define additional custom units as inputs
%   these new units may not use the name/symbol of an existing unit
%   custom units may reference built-in units, or earlier custom units
%   custom units may be combined into 1xn structure arrays
%
%   example:
%   define the density of mercury as a unit, then use this for pressures
%   Y.name = 'mecury'; Y.symb = 'hg'; Y.size = 13595.1; Y.expr = 'kgf./m3';
%   unit('in.hg','Pa') %this will fail, because "hg" is unknown
%   unit('in.hg','Pa',Y) %this will work
%
%   example:
%   define a second unit "systolic blood pressure" as 120 mmg
%   this secod unit references the first unit we added, hg.
%   Y.name = 'mecury'; Y.symb = 'hg'; Y.size = 13595.1; Y.expr = 'kgf./m3';
%   Z.name = 'systolic'; Z.symb = ''; Z.size = 120; Z.expr = 'mm.hg';
%   unit('systolic','Pa',Z,Y); %this will fail, because Z references Y
%   unit('systolic','Pa',Y,Z); %this will work
%
%   example:
%   as above, but add all new units into a single structure
%   Y.name = 'mecury'; Y.symb = 'hg'; Y.size = 13595.1; Y.expr = 'kgf./m3';
%   Z.name = 'systolic'; Z.symb = ''; Z.size = 120; Z.expr = 'mm.hg';
%   N = [Y Z];
%   unit('systolic','Pa',N); %as above, but combines custom units into one
%
%OUTPUTS
%   first output is the conversion factor from first to second input:
%   factor = unit('lb','kg'); %gives factor = 0.4536
%   if an amount is received, first output is conversion of that many:
%   amt = unit('lb','kg',[3 4 5]); %gives amt = [1.3608 1.8144 2.2680];
%
%   second output is additional information
%   [~,resu] = unit('W','kJ/hr',2)
%   resu.conv = 3.60 conversion factor
%   resu.in.orig: 'W'          %given name at input
%   resu.in.kind: 'power1.'    %kind of unit: power^1
%   resu.in.dims: 'length2.mass1.time-3.' %dimensions of unit
%   resu.in.size: 1            %input (W) divided by base unit (kgm2/s3)
%   resu.in.amt : 2            %number of units input (1 by default)
%   resu.out.orig: 'kJ/hr'     %symbol from input
%   resu.out.expr: 'energy1.time-1.'   %kind of unit: energy per time
%   resu.out.dims: 'length2.mass1.time-3.' %dimensions must match input above
%   resu.out.size: 0.2778      %input (kJ/hr) divided by base unit (kgm2/s3)
%   resu.out.amt: 7.200 %conversion factor multiplied by amount received
%
%   if custom units are received, the following fields will be added:
%   resu.customUnits.name
%   resu.customUnits.symb
%   resu.customUnits.expr
%   resu.customUnits.size
%
%DISPLAY
%   by default, numbers are displayed using the standard %g format
%   Numbers close to a multiple of pi will be shown as pi x %g
%   'eng' option display all numbers using engineering format
%   examples:
%   UNIT_CONVERT('rev/s','rad/s')          will display 2pi
%   UNIT_CONVERT('rev/s','rad/s','eng')    will display 6.28E0
%
%VERBOSITY
%   there are three levels of verbosity
%   'show'  displays the input/output unit dimesions and conversion factor
%           on by default when the function has no output(s)
%           off by default when the function has output(s)
%   'verb'  above plus shows each step of the lookup from input unis to dimension
%   'detail' also shows how each term is parsed (e.g. mm is prefix m with
%   unit m)
%
%TEMPERATURES
%   by default, temperature conversion is by scale, not absolute
%   'abstemp' option converts absolute temperatures
%   The 'abstemp' option can only be used with temperature units
%   examples:
%   unit('C','F',3)              returns 5.4 farenheit per 3 celsius
%   unit('C','F',3,'abstemp')    returns 37.4 farenheit = 3 celsius
%   unit('C/s','F/s','abstemp'); returns an error
%   The 'reltemp' option supresses the warning, for example
%   unit('C','F',3,'reltemp')    returns 5.4 with no warning
%
%MORE HELP
%unit report
%help unit>units
%help unit>outputs
%
%version: 2017.08.08
%2018-12-31 Taraz Buck - Removed "M" for mile, added "mol" for mole and "M" for molar.
%2019-03-02 Taraz Buck - Added "angstrom" for angstrom.
%% SETUP
opt = sub_setup(nargout,varargin{:}); %parse options
%% SHOW UNIT REPORT, STOP IF NO CONVERSION REQUESTED
if opt.report; sub_report(opt); end
if opt.list;   sub_list(opt);   end
if opt.help;   sub_help(opt);   end
if isequal(opt.in,'none'); return; end
%% DIMENSION LOOKUP AND CONVERSION FACTOR
%find the units in and out, and their conversion factor from base
%until all are kind.name
if opt.verb; disp('== INPUT UNIT PARSING ==='); end
opt.in = sub_lookup2kind(opt,opt.in);
opt.in = sub_lookup2dims(opt,opt.in);
if opt.verb; disp('== OUTPUT UNIT PARSING =='); end
opt.ud = sub_lookup2kind(opt,opt.ud);
opt.ud = sub_lookup2dims(opt,opt.ud);
if opt.verb; disp('============'); end
opt.conv = opt.in.size/opt.ud.size;
%check unit dimensions match
sub_check_units(opt);
%% SPECIAL HANDLING OF TEMPERATURE
%warn if converting temperature without the abstemp option
[nin,pin] = sub_str2np(opt.in.dims);
if isequal(nin,{'temperature'}) && isequal(pin,1) && ~opt.abstemp && ~opt.reltemp
    warning([mfilename ':abstemp'],'This will convert the scale of temperatures, not absolute temperatures. You may want the ''abstemp'' option.');
end
%convert absolute temperatures, if requested
if opt.abstemp; opt = sub_temp_convert(opt); end
%% ASSIGN OUTPUT
if nargout
    resu.out = rmfield(opt.ud,{'nam','pow'});
    resu.in  = rmfield(opt.in,{'nam','pow'});
    resu.conv = opt.conv;
    if opt.abstemp
        resu.in.amt = opt.in.size;
        resu.out.amt = opt.ud.size;
    else
        resu.in.amt = opt.amt;
        resu.out.amt = opt.amt.*resu.conv;
    end
    if ~isempty(opt.newunit); resu.customUnits = opt.newunit; end
    varargout{1} = resu.out.amt;
    varargout{2} = resu;
    varargout{3} = orderfields(opt);
end
if ~nargout || opt.show; sub_show_units(opt); end
end
%% SUBROUTINES : SETUP
function opt = sub_setup(nout,varargin)
%extract options from input
opts = {'show','verb','detail','report','list','help','abstemp','reltemp','eng'};
for i=1:numel(opts)
    opt.(opts{i}) = false;
end
%by default, show is true if no outputs from unit
if ~nout; opt.show = true; end
opt.newunit = [];
keep = true(size(varargin));
for i=1:numel(varargin)
    %deal with new units (defined as structures)
    if isstruct(varargin{i})
        N = varargin{i};
        rqfields = {'name','size','expr'};
        found = isfield(N,rqfields);
        if ~all(found)
            disp(rqfields(~found));
            error('new units must have the above fields');
        end
        if ~isfield(N,'symb'); N.symb = ''; end
        if isempty(opt.newunit); opt.newunit = N;
        else                     opt.newunit = [opt.newunit N];
        end
        keep(i) = false;
        continue;
    end
    
    if isnumeric(varargin{i}); continue; end %skip numbers
    switch varargin{i}
        case opts; opt.(varargin{i}) = true; keep(i) = false;
    end
end
if opt.detail; opt.verb = true; end
varargin = varargin(keep);
%any numeric inputs assigned to amt, and removed
tf = cellfun(@isnumeric,varargin);
switch sum(tf)
    case 0; opt.amt = 1;
    case 1; opt.amt = varargin{tf}; varargin = varargin(~tf);
    otherwise; error('can only accept one numeric input');
end
%remainder should be two units (strings)
switch numel(varargin)
    case 0 %only allow this if using the "report" or "list" options
        if opt.report || opt.list || opt.help; opt.in = 'none'; opt.ud = 'none';
        else error('input must contain at least one input');
        end
    case 1
        error('only one unit input received, require two');
    case 2
        opt.in = varargin{1};
        opt.ud = varargin{2};
    otherwise
        disp('BAD INPUTS');
        disp(char(varargin));
        error('can only accept up to two units');
end
%check incompatible options
if opt.abstemp && opt.reltemp; error('cannot use options ''abstemp'' and ''reltemp'' together'); end
%check inputs are correct form
if ~ischar(opt.in); error('first input must be a string'); end
if ~ischar(opt.ud); error('second input must be a string'); end
if ~isnumeric(opt.amt); error('third input must be a number or option'); end
if numel(size(opt.in))>2 || size(opt.in,1)>1 || any(size(opt.in)==0)
    error('first input must be a 1xn string, cannot be empty');
end
if numel(size(opt.ud))>2 || size(opt.ud,1)>1 || any(size(opt.ud)==0)
    error('second input must be a 1xn string, cannot be empty');
end
opt = sub_setup_prefix(opt);      %list of prefixes, e.g. M = mega = 1E6
opt = sub_setup_constant(opt);    %constants, e.g. g=9.81m/s
opt = sub_setup_kind(opt);        %kinds of unit, e.g. mass.
opt = sub_setup_unit(opt);        %all base units with scale, eg lb, kg, g
opt = sub_setup_newunit(opt);     %add custom units provided by user
if opt.report
    opt.verb = false;
    opt.detail = false;
end
end
function opt = sub_setup_kind(opt)
%creates data on "kinds" of unit
%opt = sub_setup_kind(opt)
opt.kind = sub_data_kind();
%opt.kind.dims = sub_simplify(opt.kind.expr);
%opt.kind.size = ones(size(opt.kind.dims));
%the subset of kind that is dims
opt.dims.name = opt.kind.name(opt.kind.isdim);
opt.dims.unit = opt.kind.unit(opt.kind.isdim);
opt.dims.numberof = numel(opt.dims.name);
    function kind = sub_data_kind()
        b = {...
            'mass',              'mass',              'kg',...
            'length',            'length',            'meter',...
            'time',              'time',              'second',...
            'temperature',       'temp',              'Celsius',...
            'charge',            'charge',            'coulomb',...
            'infosize',          'info',              'bit',...
            'angle',             'angle',             'rev',...
            'amount',            'amount',            'mole',...
            'area',              'length2',           'meter^2',...
            'volume',            'length3',           'meter^3',...
            'velocity',          'length/time',       'meter/s',...
            'acceleration',      'velocity/time',     'meter/second^2',...
            'force',             'mass.acceleration', 'Newton',...
            'pressure',          'force/area',        'Pascal',...
            'energy'             'force.length',      'Joule',...
            'power',             'energy/time',       'Watt',...
            'frequency',         'angle/time',        'Hertz',...
            'current',           'charge/time',       'Ampere',...
            };
        
        %removed, now defined from power&current
        %'resistance',   'R','R',    'R',...
        
        kind.name = b(1:3:end);
        kind.expr = b(2:3:end);
        kind.unit = b(3:3:end);
        kind.numberof = numel(kind.name);
        kind.isdim = false(size(kind.name));
        kind.isdim(1:8) = true;
        
        %name is a human-friendly name
        %symb is the one-letter symbol for this kind (must be unique)
        %expr is the expression for this, using symbols from "symb" field
        %unit all units are scaled relative to this, so e.g. all masses are
        %relative to kg (NOT grammes!)
        %you cannot reference this table directly when calling the
        %function, because "length/time" is not a unit, but "m/s" is.
    end
end
function opt = sub_setup_prefix(opt)
%sets up data on prefixes
%opt = sub_setup_prefix(opt)
opt.prefix.symb = 'numcHkMGTPEZY';
opt.prefix.powr = [-9 -6 -3 -2 2 3 6 9 12 15 18 21 24];
opt.prefix.name = {'nano','micro','mili','centi','hecto','kilo','mega','giga','tera','peta','exa','zetta','yotta'};
end
function opt = sub_setup_constant(opt)
%sets up conversion constants
%opt = sub_setup_con(opt)
%standard gravity, as used by e.g. lbf->N conversion
opt.con.g = 9.80665;
opt.con.abszero = -273.15;
opt.con.c = 299792458; %speed of light in m/s
opt.con.yrJ = 365.25; %julian/astronomy year
opt.con.yrcal = 365.2425; %gregorian/calendar year
%temperature conversions
opt.con.C2K = @(t) t-opt.con.abszero;
opt.con.C2R = @(t) (t-opt.con.abszero)*1.8;
opt.con.C2F = @(t) t*1.8 + 32;
opt.con.K2C = @(t) t+opt.con.abszero;
opt.con.R2C = @(t) t/1.8+opt.con.abszero;
opt.con.F2C = @(t) (t-32)/1.8;
end
function opt = sub_setup_unit(opt)
%set up data on base units
%opt = sub_setup_base(opt)
opt.unit = sub_unit_data();
opt.unit.numberof = numel(opt.unit.name);
    function out = sub_unit_data()
        b = {'meter',   'm',    1,      'length', ...
            'second',   'sec',  1,      'time', ...
            'sec',      's',    1,      'second', ...
            'celsius',  'C',    1,      'temperature', ...
            'farenheit','F',    5/9,    'celsius', ...
            'kelvin',   'K',    1,      'celsius', ...
            'rankine',  'Ra',   5/9,    'celsius', ...
            'coulomb',  'c',    1,      'charge', ...
            'amp',      'A',    1,      'coulomb/second', ...
            'ampere',    '',    1,      'amp',...
            'ohm',      'R',    1,      'watt/amp2', ...
            'volt',     'V',    1,      'watt/amp', ...
            'bit',      'bit',  1,      'infosize', ...
            'gram',     'g',    1E-3,   'mass', ...
            'ounce',    'oz',   28.3495, 'gram',...
            'revolution','rev',  1,      'angle',...
            'minute',   'min',  60,     'second',...
            'hour',     'hr',   3600,   'second', ...
            'day',      'day',  24,     'hour',...
            'year',     'yr',   opt.con.yrcal,'day',...
            'yearGregorian', '',opt.con.yrcal,'day',...
            'yearJulian','',opt.con.yrJ,'day',...
            'ton',      'T',    1E3,    'kg', ...
            'newton',   'N',    1,      'force', ...
            'litre',    '',    1E-3,   'm3', ...
            'liter',    'L',     1E-3,   'm3', ...
            'cubicCentimeter','cc',   1E-6,   'm3', ...
            'usGallon', 'gal', 1/264.172,'m3',...
            'joule',    'J',    1,      'energy', ...
            'calorie', 'cal',  4.184,  'joule', ... %lowercase cal, thermochemical calorie, at 0 degrees C. aka gram calorie.
            'foodcalorie',  'Cal'     ,4184,   'joule', ... %uppercase Cal, food calorie (aka kilocal)
            'watt',     'W',    1,      'power', ...
            'horsePowerImp','hp(I)',550,'lbf.foot/s',... %imperial/mechanical horsepower
            'horsePowerMet','hp(m)',75, 'kgf.m/s',...
            'horsepower','hp',  75,     'kgf.m/s',... %DIN 66036 (metric) version
            'herz',     'Hz',   1,      'frequency', ...
            'rpm',      '',     1,      'rev/min',...
            'degree',   'deg',  1/360,  'rev',...
            'radian',   'rad',  1/2/pi, 'rev',...
            'minuteArc',''      ,1/60,  'deg',...
            'secondArc',''      ,1/60,  'minuteArc',...
            'pascal',   'Pa',   1,      'pressure', ...
            'pound',    'lb',   0.45359237,'kg', ...
            'slug', '',         1,      'lbf.s2/ft',...
            'kiloForce', 'kgf', opt.con.g, 'N',...
            'poundForce','lbf', 0.45359237*opt.con.g,'N', ...
            'ounceForce','ozf', 28.3495E-3*opt.con.g,'N',...
            'foot',     'ft',   0.3048, 'meter',...
            'feet',     '',     0.3048, 'meter',...
            'inch',     'in',   0.0254, 'meter',...
            'yard',     '',     0.9144, 'meter',...
            'wattHour','Wh',   3600,   'joule',...
            'mile',     '',   1609.344,'meter',...
            'nauticalMile','NM',1852   ,'meter',...
            'stoke',    'St',   0.01,    'meter2/sec',...
            'psi',      'psi',  6894.75729, 'Pa',...
            'stdatmo',  'atm',  101325,'Pa',...
            'inchesMercury','inHg',3386,'Pa',...
            'bar',      'bar',1E5,'Pa',...
            'gee',      '',opt.con.g,'acceleration',...
            '',         'kph',1/3.6,'m/s',...
            'lightSpeed','cee',opt.con.c,'velocity',...
            'lightYear','LY',1,'lightSpeed*yearJulian',...
            'acre',     '',66*660,'ft2',...
            'hectare',  '',1E4,'m2',...
            'henry',    'H',1,'weber/ampere',...
            'siemen','',1,'ohm-1',...
            'weber','Wb',1,'volt/second',...
            'tesla','',1,'weber/area',...
            'gauss','G',1E-4,'tesla',...
            'mole',     'mol',      1,                'amount', ...
            'Molar',    'M',        1,                'mole./liter', ...
            'molecule', 'molecule', 1/6.02214076E23,  'mole', ...
            'angstrom', '',         1E-10,            'meter',...
            };
        out.name = b(1:4:end);
        out.symb = b(2:4:end);
        out.size = [b{3:4:end}];
        out.expr = b(4:4:end);
    end
end
function opt = sub_setup_newunit(opt)
%adds new (custom) units to the set of units
U = opt.unit;
pat = 'custom unit "%s" bad format: field "%s" must be a %s\n';
pat2 = 'custom unit "%s" clashes: field "%s" value "%s" already exists\n';
for ii = 1:numel(opt.newunit)
    new = opt.newunit(ii);
    if ~isfield(new,'symb'); new.symb = ''; end
    %check the new unit is compatible
    if ~ischar(new.name);       error(pat,'','name','string'); end
    if ~ischar(new.expr);       error(pat,new.name,'expr','string'); end
    if ~ischar(new.symb);       error(pat,new.name,'symb','string'); end
    if ~isnumeric(new.size);    error(pat,new.name,'size','number'); end
    
    %warn about anything being overwritten
    if ~isempty(new.name) && ismember(new.name,U.name); error(pat2,new.name,'name',new.name); end
    if ~isempty(new.symb) && ismember(new.symb,U.symb); error(pat2,new.name,'symb',new.symb); end
    
    %ud = sub_base_unit(opt,new.expr);
    U.name{end+1} = new.name;
    U.symb{end+1} = new.symb;
    U.size(end+1) = new.size;
    U.expr{end+1} = new.expr;
end
opt.unit = U;
if (numel(opt.newunit) && opt.show) || opt.verb
    fprintf('Using %d custom units:\n',numel(opt.newunit));
    pat = 'name "%s" expression 1 %s = %g x %s\n';
    for ii=1:numel(opt.newunit)
        fprintf(pat,...
            opt.newunit(ii).name,...
            opt.newunit(ii).symb,...
            opt.newunit(ii).size,...
            opt.newunit(ii).expr);
    end
    fprintf('\n');
end
end
%% PARSING OF UNITS AND DIMENSION LOOKUP
function in = sub_lookup2dims(opt,in)
if opt.verb; disp('--- start lookup from kind to dimensions ---'); end
if ~isfield(in,'pow')
    [in.nam,in.pow] = sub_str2np(in.kind);
end
if ~isfield(in,'let')
    in.let = sub_term_type(opt,in.kind);
end
v.let = in.let;
v.pow = in.pow;
v.nam = in.nam;
v.siz = nan;
in = rmfield(in,{'let','pow','nam'});
while any(v.let~='A')
    if ~all(v.let=='A' | v.let=='B'); keyboard; end
    I = find(v.let=='B',1,'first');
    [tf,loc] = ismember(v.nam{I},opt.kind.name);
    if ~tf; keyboard; end
    str = opt.kind.expr{loc};
    if opt.verb
        fprintf(['replaced kind with expression "' v.nam{I} '" -> "' str '"\n']);
    end
    [n,p] = sub_str2np(str);
    let = sub_term_type(opt,str);
    v.nam = [v.nam n];
    v.pow = [v.pow p*v.pow(I)];
    v.let = [v.let let];
    v.pow(I) = 0;
    v = sub_simple(v);
    if opt.verb; fprintf('%s\n',sub_np2str(v.nam,v.pow)); end
end
if opt.verb; disp('--- lookup to dimensions completed ---'); end
[v.nam,order] = sort(v.nam);
v.pow = v.pow(order);
in.dims = sub_np2str(v.nam,v.pow);
[tf,loc] = ismember(v.nam,opt.kind.name);
if ~all(tf); keyboard; end
in.nam = v.nam;
in.pow = v.pow;
in.unit = sub_np2pretty(opt.kind.unit(loc),v.pow);
end
function in = sub_lookup2kind(opt,str)
in.orig = str;
v.siz = 1;
[v.nam,v.pow] = sub_str2np(str);
v.let = sub_term_type(opt,str);
if opt.verb
    fprintf('Starting unit %s\n',in.orig);
end
while ~all(v.let=='A' | v.let=='B')
    if any(v.let=='I'); keyboard; end
    v = sub_rem_prefix(opt,v);
    v = sub_symb2name(opt,v);
    v = sub_name2expr(opt,v);
    v = sub_simple(v);
end
if opt.verb; disp('--- lookup to kind completed ---'); end
[v.nam,order] = sort(v.nam);
v.pow = v.pow(order);
v.let = v.let(order);
in.kind = sub_np2str(v.nam,v.pow);
in.size = v.siz;
in.nam = v.nam;
in.pow = v.pow;
in.let = v.let;
end
function v = sub_name2expr(opt,v)
while any(v.let=='C')
    I = find(v.let=='C',1,'first');
    [tf,loc] = ismember(v.nam{I},opt.unit.name);
    if ~tf; keyboard; end
    scale = opt.unit.size(loc)^v.pow(I);
    v.siz = v.siz*scale;
    %split the new expression into terms and append
    if opt.verb; old = v.nam{I}; end
    str = opt.unit.expr{loc};
    [n,p] = sub_str2np(str);
    let = sub_term_type(opt,str);
    v.nam = [v.nam n];
    v.pow = [v.pow p*v.pow(I)];
    v.let = [v.let let];
    %remove the term we just replaced
    v.nam(I) = [];
    v.pow(I) = [];
    v.let(I) = [];
    v = sub_simple(v);
    if opt.verb
        fprintf('replaced %s with expression %s\n',old,str);
        fprintf('%s x %e\n',sub_np2str(v.nam,v.pow),v.siz);
    end
end
end
function out = sub_simple(v)
keep = v.pow~=0;
v.let = v.let(keep);
v.pow = v.pow(keep);
v.nam = v.nam(keep);
[out.nam,I,J] = unique(v.nam);
for i=1:numel(out.nam)
    out.pow(i) = sum(v.pow(J==i));
end
if isempty(out.nam)
    out.pow = [];
    out.let = '';
end
out.let = v.let(I);
    keep = out.pow~=0;
    out.let = out.let(keep);
    out.pow = out.pow(keep);
    out.nam = out.nam(keep);
    out.siz = v.siz;
end
function v = sub_rem_prefix(opt,v)
%strip full-name prefixes
for I = find(ismember(v.let,'EF'))
    [prestr,pow10] = sub_find_prefix(opt,v.nam{I},'name');
    v.nam{I} = regexprep(v.nam{I},['^' prestr],'');
    scale = 10^(pow10*v.pow(I));
    v.siz = v.siz * scale;
    v.let(I) = sub_term_type(opt,v.nam{I});
    v = sub_simple(v);
end
%strip letter prefixes
for I = find(ismember(v.let,'GH'))
    [prestr,pow10] = sub_find_prefix(opt,v.nam{I},'symb');
    if opt.verb; fprintf('removing prefix %s\n',prestr); end
    v.nam{I} = regexprep(v.nam{I},['^' prestr],'');
    scale = 10^(pow10*v.pow(I));
    v.siz = v.siz * scale;
    v.let(I) = sub_term_type(opt,v.nam{I});
    v = sub_simple(v);
end
if opt.verb; fprintf('%s x %g\n',sub_np2str(v.nam,v.pow),v.siz); end
    function [strud,pow10] = sub_find_prefix(opt,str,name_or_symb)
    switch name_or_symb
        case 'name'; from = opt.prefix.name;
        case 'symb'; from = cellstr(opt.prefix.symb');
    end
    to = opt.prefix.powr;
    for i=1:numel(from)
        pat = ['^' from{i}];
        if ~isempty(regexp(str,pat,'once'))
            strud = from{i};
            pow10 = to(i);
            return
        end
    end
    if ~exist('strud','var'); keyboard; end
    end
end
function v = sub_symb2name(opt,v)
while any(v.let=='D')
    if any(v.let=='I'); keyboard; end
    I = find(v.let=='D');
    [tf,loc] = ismember(v.nam(I),opt.unit.symb);
    if opt.verb
        for i=1:numel(I)
            fprintf(...
                'replacing symb %s with name %s\n',...
                v.nam{I(i)},opt.unit.name{loc(i)});
        end
    end
    if ~all(tf); keyboard; end
    v.nam(I) = opt.unit.name(loc);
    v.let(I) = 'C';
    v = sub_simple(v);
    if opt.verb; fprintf('%s x %e\n',sub_np2str(v.nam,v.pow),v.siz); end
end
end
function let = sub_term_type(opt,str)
%let = sub_term_type(opt,str)
%assign each term a letter code determining its type
%if more than one code applies, the first is used
%A dims.name
%B kind.name
%C base.name
%D base.symb
%E prefix name + base.name
%F prefix name + base.symb
%G prefix symb + base.name
%H prefix symb + base.symb
%I other (error)
if iscell(str); n = str;
else;           n = sub_str2np(str);
end
let = repmat('I',[1 numel(n)]);
for i=1:numel(n)
    let(i) = sub_term_single(opt,n{i});
end
if isempty(n); let = ''; return; end
if any(let=='I')
    disp('BAD UNITS');
    fprintf(1,'"%s"\n',n{let=='I'});
    error('could not parse the units shown above');
end
    function letter = sub_term_single(opt,str)
        %catches in order:
        %dims.name, kinds.name, base.name, base.symb
        if ismember(str,opt.kind.name(opt.kind.isdim)); letter = 'A';
            if opt.detail; fprintf('%s is kind %s : a base dimension\n',str,letter); end
            return;
        end
        if ismember(str,opt.kind.name); letter = 'B';
            if opt.detail; fprintf('%s is kind %s : a kind not base dimension\n',str,letter); end
            return;
        end
        if ismember(str,opt.unit.name); letter = 'C';
            if opt.detail; fprintf('%s is kind %s : a unit name\n',str,letter); end
            return;
        end
        if ismember(str,opt.unit.symb); letter = 'D';
            if opt.detail; fprintf('%s is kind %s: a unit symbol\n',str,letter); end
            return;
        end
        
        %catches prefix.symb + name and prefix.symb + symb
        %create pattern to find all scale prefixes
        pat = sub_pat_prefix(opt.prefix.name);
        pre = regexp(str,pat,'match');
        if ~isempty(pre)
            for j=1:numel(pre) %for each possible prefix
                str = regexprep(str,['^' pre{j}],'');
                if isempty(str)
                    fprintf('%s is a prefix with no unit, can''t use this\n',[pre{j} str]);
                    continue
                end
                if ismember(str,opt.unit.name); letter = 'E';
                    if opt.detail; fprintf('%s is kind %s : prefix %s + name %s\n',[pre{j} str],letter,pre{j},str); end
                    return;
                end
                if ismember(str,opt.unit.symb); letter = 'F';
                    if opt.detail; fprintf('%s is kind %s : prefix %s + symbol %s\n',[pre{j} str],letter,pre{j},str); end
                    return;
                end
            end
        end
        
        pat = sub_pat_prefix(opt.prefix.symb);
        pre = regexp(str,pat,'match');
        
        if ~isempty(pre)
            for j=1:numel(pre) %for each possible prefix
                str = regexprep(str,['^' pre{j}],'');
                if isempty(str) 
                    fprintf('%s is a prefix with no unit, can''t use this\n',str);
                    continue;
                end
                if ismember(str,opt.unit.name); letter = 'G';
                    if opt.detail; fprintf('%s is kind %s : prefix %s + symbol %s\n',[pre{j} str],letter,pre{j},str); end
                    return;
                end
                if ismember(str,opt.unit.symb); letter = 'H';
                    if opt.detail; fprintf('%s is kind %s : prefix %s + symbol %s\n',[pre{j} str],letter,pre{j},str); end
                    return;
                end
            end
        end        
        letter = 'I'; return
    end
    function pat = sub_pat_prefix(cstr)
        if ~iscell(cstr); cstr = cellstr(cstr'); end
        pat = ['^(' sprintf('%s|',cstr{:}) ')'];
    end
end
function [n,pow] = sub_str2np(str)
%splits string str into units n and powers pow
%[n,pow] = sub_str2np(str)
%split into terms
%each term may start with / . * or whitespace \s
%then a number of letters
%then some other stuff until the next match
%TODO get non-greedy flag to work with regexp, so get all terms with one
%call.
if isempty(str); n = {}; pow = []; return; end
%if string starts with 1/unit, replace with /unit
pat = '^1/[a-z|A-Z|\(|\)]'; if regexp(str,pat); str(1) = ''; end
pat = '[\.|\*|/|\s+]?[a-z|A-Z|\(|\)]+';
beg = regexp(str,pat,'start');
if beg(1)~=1; error('unit should start with a letter or /'); end
bak = unique([beg(2:end)-1 numel(str)]);
for i=1:numel(beg)
    terms{i} = str(beg(i):bak(i)); %#ok<AGROW>
end
%catch use of - as separator, discard the -
terms = regexprep(terms,'-$','');
%TODO here : add parsing of bracketed terms, so
%kg/m/s2   = kg1m-1s-2
%kg/(m/s2) = kg1m-1s+2
%if a term starts with /, invert the power of all following terms
%e.g. kg/m.s is kg1m-1s-1.
%CHANGE 2016.03.18 if a term start with / invert the power of just the
%following term, not all following terms. so example above
%kg/m.s is now kg1m-1s1
pat = '^/';
pow = double(cellfun('isempty',regexp(terms,pat)));
pow(pow==0) = -1;
%pow = cumprod(pow);    %this changes the "all following terms" part.
terms = regexprep(terms,pat,'');
%remove leading . and * from each term
pat = '^[\.|\*|\s+]';
terms = regexprep(terms,pat,'');
%parse and remove the units in each term
pat = '^[a-z|A-Z|\(|\)]+';
n = regexp(terms,pat,'match','once');
terms = regexprep(terms,pat,'');
%each term should now just be a power.
%remove starting ^ before powers
terms = regexprep(terms,'^\^','');
%remove ending .
terms = regexprep(terms,'\.$','');
%parse powers term-by-term to catch fractions (/)
for i=1:numel(terms)
    if isempty(terms{i}); val(i) = 1; continue; end %#ok<AGROW>
    if regexp(terms{i},'/')
        a = regexp(terms{i},'(?<num>.+)/(?<den>.+)','names');
        val(i) = str2double(a.num)/str2double(a.den); %#ok<AGROW>
    else
        val(i) = str2double(terms{i}); %#ok<AGROW>
    end
end
pow = pow .* val;
%consolidate multiples of the same unit
if numel(unique(n))<numel(n)
    
    [u,~,c] = unique(n);
    p = [];
    for i=1:numel(u);
        p(i) = sum(pow(c==i)); %#ok<AGROW>
    end
    n = u;
    pow = p;
end
if any(isnan(pow))
    I = find(isnan(pow),1,'first');
    try
        error('failed to parse %s, cannot parse power of term %s',str,n{I})
    catch
        keyboard
    end
    keyboard;
end
end
function str = sub_np2str(n,p)
str = '';
for i=1:numel(n)
    str = [str sprintf('%s%g.',n{i},p(i))];
end
end
function str = sub_np2pretty(n,p)
num = sub_np2nice(n(p>0),p(p>0)); if isempty(num); num = 1; end
str = num;
den = sub_np2nice(n(p<0),-p(p<0)); if isempty(den); return; end
str = [str '/(' den ')'];
    function s = sub_np2nice(a,b)
    s = '';
    for i = 1:numel(a)
        try
            s = [s ' x ' a{i}];
        catch
            keyboard
        end
        if b(i)==1; continue; end
        s = [s '^' num2str(b(i))];
    end
    if numel(s); s(1:3) = ''; end
    end
end
function strud = sub_simplify(str,flag_verb)
%simplifies the units of str, removing zero-power units
%str = sub_simplify(str,flag_verb)
if nargin<2; flag_verb = false; end
%if cell array, recursively call self
if iscellstr(str)
    for i=1:numel(str)
        strud{i} = sub_simplify(str{i},flag_verb); %#ok<AGROW>
    end
    return
end
%simplify expressions
[n,pow] = sub_str2np(str);
%remove any with zero power and consolidate into unique names
keep = pow~=0;
pow = pow(keep);
n = n(keep);
%create string out
strud = '';
if ~isempty(n)
    for u = unique(n)
        [~,loc] = ismember(u,n);
        p = sum(pow(loc));
        strud = [strud u{1} num2str(p) '.']; %#ok<AGROW>
    end
end
if flag_verb; fprintf(1,'simplify "%s" -> "%s"\n',str,strud); end
end
function opt = sub_temp_convert(opt)
%check input is purely temperature
%unit check has been performed by this point, so know input/output match
%check unit is just temperature
pat = '^temperature[\d|\.]+$';
if isempty(regexp(opt.in.dims,pat,'once')); error('units must be just temperatures when using abstemp option'); end
if isempty(regexp(opt.ud.dims,pat,'once')); error('units must be just temperatures when using abstemp option'); end
%check unit is a single temperature (no C/F C*K and so on)
[nin,pin] = sub_str2np(opt.in.orig);
[nud,pud] = sub_str2np(opt.ud.orig); %#ok<NASGU>
if ~isequal(pin,1); error('can only use the "abstemp" option with temperatures to the power of 1'); end
if numel(nin)>1; error('absolute temperature conversion can only accept a single unit at a time'); end
if numel(nud)>1; error('absolute temperature conversion can only accept a single unit at a time'); end
nin = nin{1}; nud = nud{1};
%express all in Celsius
switch lower(nin);
    case {'c','celsius'};   tin = opt.amt;
    case {'k','kelvin'};    tin = opt.con.K2C(opt.amt);
    case {'f','farenheit'}; tin = opt.con.F2C(opt.amt);
    case {'ra','rankine'};  tin = opt.con.R2C(opt.amt);
    otherwise; error('unknown temperature unit');
end
%convert to output units
switch lower(nud)
    case {'c','celsius'};   tud = tin;
    case {'k','kelvin'};    tud = opt.con.C2K(tin);
    case {'f','farenheit'}; tud = opt.con.C2F(tin);
    case {'ra','rankine'};  tud = opt.con.C2R(tin);
    otherwise; error('unknown temperature unit');
end
opt.in.size = opt.amt;
opt.ud.size = tud;
end
%% OUTPUT
function sub_check_units(opt)
%check that input/output units are the same
[nin,pin] = sub_str2np(opt.in.dims);
[nud,pud] = sub_str2np(opt.ud.dims);
if ~isequal(nin,nud)
    sub_disp_fromto(opt);
    error('dimensions in/out do not match');
end
if isempty(nin); return; end
%warn or error on power mismatch
err = abs(pin - pud);
worst = max(err);
tol = [1E-6 1E-3];
%if all within 1E-6, don't care
if worst<tol(1); return; end
%if any worse than 1E-3, error
if worst>tol(2)
    sub_disp_fromto(opt);
    error(['unit powers do not match to within ' num2str(tol(2)) ', worst error ' num2str(worst)]);
end
%warn about slight mismatch, but continue
if worst>=tol(1) && worst<=tol(2)
    sub_disp_fromto(opt);
    warning(['unit powers do not match, worst error ' num2str(worst)])
end
    function sub_disp_fromto(opt)
        disp(['from unit: ' opt.in.orig])
        disp(['  to unit: ' opt.ud.orig]);
        disp(' ');
        disp(['from kind: ' opt.in.kind]);
        disp(['  to kind: ' opt.ud.kind]);
        disp(' ');
        disp(['from dimensions: ' opt.in.dims]);
        disp(['  to dimensions: ' opt.ud.dims]);
    end
end
function str = sub_pretty_number(opt,num)
%engineering format, if requested
%to nearest 1E3 below
if opt.eng
    pow = floor(log10(num)/3)*3;
    str = sprintf('%ge%d',num/(10^pow),pow);
    return
end
%is it a round fraction ?
[n,d] = rat(num);
if d==1 && n>0; str = sprintf('%d',n);  return; end
if n==1;        str = sprintf('%d/%d',n,d); return; end
%is it a pi fraction ?
[nr,dr] = rat(num/pi);
if dr==1 && nr==1;      str = 'pi';                     return; end
if dr==1 && nr>1;       str = sprintf('%gpi',nr);       return; end
if nr==1;               str = sprintf('pi/%d',dr);      return; end
%is it a pi multiple ?
[nr,dr] = rat(num*pi);
if nr==1 && dr==1; str = '1/pi';         return; end
if dr==1 && nr>1 ; str = sprintf('%d/pi',nr);   return; end
if nr==1;   str = sprintf('1/%dpi',dr); return; end
%default, g format (if none of the others applied
str = sprintf('%G',num);
end
function sub_show_units(opt)
if opt.abstemp
    for i=1:numel(opt.in.size)
        fprintf(1,'%s %s\t=\t%s %s\n',...
            sub_pretty_number(opt,opt.in.size(i)),...
            opt.in.orig,...
            sub_pretty_number(opt,opt.ud.size(i)),...
            opt.ud.orig);
    end
    return
end
v = opt.conv.*opt.amt;
for i=1:numel(opt.amt)
    fprintf(1,'%s %s\t=\t%s %s\n',...
        sub_pretty_number(opt,opt.amt(i)),...
        opt.in.orig,...
        sub_pretty_number(opt,v(i)),...
        opt.ud.orig);
end
disp('-');
fprintf(1,'%s\t=\t%s\tx %s\n',opt.in.orig,opt.ud.orig,sub_pretty_number(opt,opt.conv));
fprintf(1,'%s\t=\t%s\tx %s\n',opt.ud.orig,opt.in.orig,sub_pretty_number(opt,1/opt.conv));
disp('-');
sub_show_dimensions(opt,opt.in,'from');
disp('-');
sub_show_dimensions(opt,opt.ud,'to');
function sub_show_dimensions(opt,in,str)
if nargin<3; str = ' '; end
[n,p] = sub_str2np(in.dims);
ns = opt.kind.name(opt.kind.isdim);
us = opt.kind.unit(opt.kind.isdim);
ss = opt.kind.name(opt.kind.isdim);
ps = zeros(size(ns));
[tf,loc] = ismember(n,ns);
if ~all(tf); keyboard; end
ps(loc) = p;
fprintf(1,'%-6s: %s\n',str,in.orig);
    fprintf(1,'scale : x%s %s\n',...
        sub_pretty_number(opt,in.size),...
        in.unit);
    if isempty(in.kind)
        fprintf(1,'%s : %s\n','kind','nondimensional');
        fprintf(1,'%s : %s\n','dims','nondimensional');
    else
        fprintf(1,'%s : %s\n','kind ',in.kind);
        fprintf(1,'%s : %s\n','dims ',in.dims);
    end
end
end
function sub_report(opt)
clc;
disp('==================')
disp('== UNIT REPORT ===');
disp('==================');
fprintf('These are the inputs recognised by the function %s\n',mfilename);
disp('================');
disp('=== OPTIONS ====');
disp('================');
disp(char(setxor(fieldnames(opt),{'amt','prefix','con','kind','dims','unit','in','ud','newunit'})));
disp(' ');
disp('================')
disp('=== PREFIXES ===');
disp('================')
for i=1:numel(opt.prefix.symb)
    fprintf('%c %7s 10E%+d\n',...
        opt.prefix.symb(i),...
        opt.prefix.name{i},...
        opt.prefix.powr(i));
end
disp('Each term in a unit expression is allowed a single prefix')
disp('For example: W, kW, MW, GW, nanoW, megaW')
disp(' ');
disp('======================');
disp('=== KINDS OF UNITS ===');
disp('======================');
fprintf('%12s %-5s %-14s\n',...
    'name','dim?','unit');
ny = {'no','yes'};
ny = ny(opt.kind.isdim+1);
for i=1:numel(opt.kind.name)
    fprintf('%12s %-5s %-14s\n',...
        opt.kind.name{i},...
        ny{i},...
        opt.kind.unit{i})
end
disp('all input units are reduced to a combination of kinds of unit, and then dimensions');
disp('these are checked against each other to check the from/to dimensions match');
disp('dim? shows if this is a basic dimension')
disp('for example "length" is a basic dimension, but "area" is not, as area has dimension length^2')
disp('all units are reduced to basic dimensions to check if dimensions match');
disp(' ');
disp('========================')
disp('=== UNITS (BY NAME) ====');
disp('========================')
fprintf('%16s %5s %15s\n',...
    'NAME','SYMB','EXPRESSION');
%sort units by name
[~,order] = sort(lower(opt.unit.name));
for i=1:numel(opt.unit.name)
    fprintf('%16s %5s %8s x %s\n',...
        opt.unit.name{order(i)},...
        opt.unit.symb{order(i)},...
        sub_pretty_number(opt,opt.unit.size(order(i))),...
        opt.unit.expr{order(i)});
end
disp(' ');
disp('=========================')
disp('=== UNITS (BY SYMBOL) ===')
disp('=========================')
fprintf('%16s %5s %15s\n',...
    'NAME','SYMB','EXPRESSION');
%sort units by name
[~,order] = sort(lower(opt.unit.symb));
for i=1:numel(opt.unit.name)
    fprintf('%16s %5s %8s x %s\n',...
        opt.unit.name{order(i)},...
        opt.unit.symb{order(i)},...
        sub_pretty_number(opt,opt.unit.size(order(i))),...
        opt.unit.expr{order(i)});
end
%opt.in = sub_lookup2kind(opt,opt.in);
%opt.in.symb = opt.in.orig;
%opt.in = sub_lookup2dims(opt,opt.in);
%20151021 reintroduce by dimesions
%find basic dimensions of each kind
for i=1:numel(opt.kind.name)
    v.kind = opt.kind.name{i};
    v = sub_lookup2dims(opt,v);
    k{i} = v.dims;
    clear('v');
end
%find basic dimensions of each base unit
for i=1:numel(opt.unit.name)
    v = sub_lookup2kind(opt,opt.unit.name{i});
    v = sub_lookup2dims(opt,v);
    u{i} = v.dims;
end
disp(' ');
disp('=======================')
disp('=== UNITS (BY TYPE) ===')
disp('=======================')
fprintf('%16s %5s   %s\n','NAME','SYMB','EXPRESSION');
used = false(size(u));
for i=1:numel(opt.kind.name)
    fprintf('--- %s ---\n',opt.kind.name{i});
    for j=1:numel(opt.unit.name)
        if isequal(k{i},u{j})
            fprintf('%16s %5s %8s x %s\n',...
                opt.unit.name{j},...
                opt.unit.symb{j},...
                sub_pretty_number(opt,opt.unit.size(j)),...
                opt.unit.expr{j});
            used(j) = true;
        end
    end
end
fprintf(1,'--- %s ---\n','other');
for j=find(~used)
    fprintf('%16s %5s %8s x %s\n',...
        opt.unit.name{j},...
        opt.unit.symb{j},...
        sub_pretty_number(opt,opt.unit.size(j)),...
        opt.unit.expr{j});
end
disp(' ');
disp('all your input units should be from the "SYMB" or "NAME" column above.');
disp(' ');
disp('terms must be separated by one of *-./ e.g. N-m or m/s')
disp('each term may have zero or one prefixes, e.g. kN/m or kg/cm3')
disp('each term may have zero or one power, a number, e.g. kg/m3 or lbf*in-2.')
disp('powers may be fractions, e.g. kg^1/2 or mm1/4');
disp(' ');
disp('The expressions in this table are mostly references to the "name" column of other units');
disp('a few are references to the kind of unit. These are the units that define that kind');
disp('for example, "meter" has an expression of "length", showing it is the defining unit for kind "length"');
end
function sub_list(opt)
a = [opt.unit.name opt.unit.symb];
tf = cellfun(@isempty,a);
a = a(~tf);
a = unique(a);
[~,order] = sort(lower(a));
a = a(order);
disp('=========================');
disp('=== LIST OF ALL UNITS ===');
disp('=========================');
disp(char(a));
disp('=========================');
end
%% HELP SUBROUTINES
%these functions don't do anything, an shouldn't be called except via help
%e.g. help unit>options
%this section is still WIP
function outputs()
%== DETAILED OUTPUT HELP ==
%the function uses nargout to determine the number of outputs requested
%
%ZERO OUTPUTS
%if no outputs are requested, nothing is returned
%instead, the conversion between the two units is shown on the command line
%
%ONE OUTPUT
%the first output is the conversion between the two units, for example:
%a = unit('ft','m',[3 4 5]);
%returns a = [0.9144    1.2192    1.5240], because
%3 ft = 0.9144 m
%4 ft = 1.2194 m
%etc
%
%if no numerical input is recievied, it is assumed you want to conver one
%of the first unit, for example:
%a = unit('ft','m');
%returns a = 0.3048, because
%1 foot = 0.3048 meters
%
%TWO OUTPUTS
%The second output contains additional information about the units
%[~,b] = unit('N','kg.m/s2',2)
%b is a structure with fields:
%b.conv     conversion factor between the two units
%b.in       details of the first (input) unit
%b.out      details of the second (output) unit
%b.customUnits details of any custom units received as an input
%
%the 'in' and 'out' fields have subfields:
%the full list of fields is:
%b.conv = 3.60 conversion factor
%b.in.orig:  'N'             %name of unit from input
%b.in.kind:  'force1.'       %kind and power of unit
%b.in.dims:  'length1.mass1.time-2.' %dimensions of unit
%b.in.size:  1               %size of unit relative to base
%b.in.amt:   2    %amount of unit to convert
%b.out.orig: 'lb*gee' %name of unit from input (pounds * gee acceleration)
%b.out.kind: 'length1.mass1.time-2.' %kind and power of unit
%b.out.dims: 'length1.mass1.time-2.' %base dimensions
%b.out.size: 4.4482 %size of unit relative to base dimension
%b.out.amt:  0.4496 %the input amount, converted to this unit
%b.customUnits.name
%b.customUnits.symb
%b.customUnits.expr
%b.customUnits.size
%
%THREE OUTPUTS
%The third output is the full details of all values used. This is useful
%for documentation and validation purposes
%[~,~,c] = unit('N',kg.m/s');
end
function units()
%=== UNIT DEFINITIONS ==
%units are divided into three tiers: unit, kind and dims
%
%UNIT is the full list of all units. There are currently 66 of these
%when used, these are reduced to combinations of the 18 KINDs of unit
%Finally, they are reduced to combinations of the 8 BASE units
%Example:
%Litre is a unit, of kind "volume" and base units "length^3"
%unit('report') shows a report of all units, kinds and dimensions
%
%Both the inputs are converted to base, then checked to see if they match.
%The UNIT and KIND need not match, but the BASE must.
%example: unit('gee','lbf/lb')
%kind does not match, they are "acceleration" and "force per mass"
%base does match: length/time2 both
%
%BASE
%There are 8 base units
%mass        (kilogram)
%length      (meters)
%time        (seconds)
%temperature (celsius)
%charge      (coulomb)
%infosize    (bit)
%rev         (revolutions) used in e.g. rpm and Hz
%amount      (moles)
%
%This sort of analysis is typically called a "dimensional analysis", and
%the base units above would be called dimensions. I have chosen to call them
%base units instead, as people can get touchy about what is and is not a 
%proper "dimension". You're free to call them what you want :)
%
%each base unit has an si unit attached. This is the reference magnitude of
%that base.
%
%KIND
%kinds of unit are the 7 dimensions above, plus 10 others.
%The 10 added are defined in terms of each other, and base units:
%area           length^2
%volume         length^3
%velocity       length/time
%acceleration   velocity/time
%force          mass.acceleration
%pressure       force/area
%energy         force.length
%power          energy/time
%frequency      angle/time
%current        charge/time
%
%kinds are helpful for writing the expression of units - it's
%more intuitive to define a pressure as force/area than mass/length/time2
%
%UNIT
%this is the full list of all 62 units currently defined in the code -
%althought you also have the option to add more at run-time.
%each unit is defined by four fields:
%name : long "proper" name of the unit
%symb : short "symbolic" name
%expr : expression for the unit that allows you to reduce it to base
%size : scale of the unit relative to expression.
%
%the "expression" field is used to reduce the unit to kinds, and eventually base units.
%for example:
%PoundForce has expression 4.44xN
%that means any instances of poundforce will be replaced with 4.44 lots of
%N (Newtons).
%from the same table, newton is a unit with expression 1xforce.
%that means any instances of newton will be replaced with the kind "force".
%
%DEFINING UNITS
%It's simple enough to add new units, either by editing the code, or by
%adding them as an input at runtime. There are a few restrictions on the
%naming:
%
%The name and symbol may only contain letters a-z, A-Z and brackets.
%No letters, whitespace, punctuation, etc. are allowed.
%
%names are typically lowercase, with compound names having capitalisation:
%e.g. newton and nauticalMile.
%
%symbols are a mix of upper-and-lowercase, preferentially lowercase.
%
%The expression "expr" must allow you to reduce it to base dimensions.
%Preferably, the expression should be a reference direct to kinds of unit,
%but sometimes it's more helpful to write it in terms of other units.
%Be really careful when doing this that you don't accidentally introduce
%circular references.
end
%% DEVNOTES
% TODO handling of dimensionless ratios, e.g. decibel
