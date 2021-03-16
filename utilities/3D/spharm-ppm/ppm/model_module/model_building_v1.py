from spp_v3 import *
# sys.argv[0]

temp_folder_name='/Users/piepi/CellOrganizer/cellorganizer3/demos/3D/demo3Dpoint_process_module/temp/'
folder_name='/Users/piepi/CellOrganizer/Josh/human brain imaging 2nd round/'
regex_input='.*40x.*.nd2'
files=os.listdir(folder_name)
regex = re.compile(regex_input)
img_names=sorted_nicely(filter(regex.search,files))
quads = []
for img_name in img_names:
    pts = np.loadtxt(temp_folder_name+img_name+'__centers.txt', dtype=float)
    maxx = 2048
    maxy = 2048
    maxz = find_zmax(temp_folder_name,img_name)
    minx = 1
    miny = 1
    minz = 1
    pts_scale = np.column_stack(((pts[:,0]-minx)/(maxx-minx),(pts[:,1]-miny)/(maxy-miny),(pts[:,2]-minz)/(maxz-minz)))

    cur_pat = PointPattern(pts_scale,dimension=3)
    print 'PointPattern complete'
    quads.append(QuadScheme(cur_pat,50))
    print 'QuadScheme complete'

betas = []
failed = []
for qind in xrange(0,len(quads)):
    q = quads[qind]
    try:
        factor_matrix = np.loadtxt(temp_folder_name+img_names[qind]+'__factor_matrix.txt')
        beta = FitModel(q,temp_folder_name,img_names[qind],factor_matrix)
        betas.append(beta)
    except:
        failed.append(qind)
        print 'failed',qind
        pass

with open('coefficient.txt','w') as f:
    for val in betas:
        for b in val:
            f.write(str(b)+' ')
        f.write('\n')
    f.write('\n')
    mean_beta = np.mean([betas[bind] for bind in xrange(0,len(betas))],axis=0)
    for b in mean_beta:
        f.write('full model: ')
        f.write(str(b)+' ')
    f.write('\n')

fit_quads = [quad for qind,quad in enumerate(quads) if qind not in failed]

cv_ll = []
for qind in xrange(0,len(img_names)):
    held = fit_quads[qind]
    mean_beta = np.mean([betas[bind] for bind in xrange(0,len(betas)) if bind != qind],axis=0)
    logl = inhomoLogLike(held.pat,mean_beta,temp_folder_name,img_names[qind])
    cv_ll.append(logl)

    with open('likelihood.txt','a') as f:
        f.write(str(logl)+'\n')

avg_ll = np.mean(cv_ll)
var_ll = np.var(cv_ll)

with open('likelihood.txt','a') as f:
        f.write('\n')
        f.write('average: '+str(avg_ll)+'\n')
        f.write('variance: '+str(var_ll)+'\n')