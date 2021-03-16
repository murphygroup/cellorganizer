from spp_tim_v1 import *

folder_name=sys.argv[1]
regex_input=sys.argv[2]
root_path=sys.argv[3]
temp_path=sys.argv[4]
temp_folder_name=root_path+temp_path
files=os.listdir(folder_name)
regex = re.compile(regex_input)
img_names=sorted_nicely(filter(regex.search,files))
quads = []
for img_name in img_names:
    pts = np.loadtxt(temp_folder_name+'others/'+img_name+'__centers.txt', dtype=float)
    maxx,maxy,maxz=find_xmax_ymax_zmax(temp_folder_name,img_name)
    # print maxx,maxy,maxz
    minx = 1
    miny = 1
    minz = 1
    pts_scale = np.column_stack(((pts[:,0]-minx)/(maxx-minx),(pts[:,1]-miny)/(maxy-miny),(pts[:,2]-minz)/(maxz-minz)))
    cur_pat = PointPattern(pts_scale,dimension=3)
    print 'PointPattern complete'
    quads.append(QuadScheme(cur_pat,50))
    print 'QuadScheme complete'

full_model,betas,failed=FitModel_all(quads,temp_folder_name,img_names)
with open('full_model.txt','a') as f:
    for b in full_model:
        f.write(str(b)+' ')
    f.write('\n')

fit_quads = [quad for qind,quad in enumerate(quads) if qind not in failed]

cv_ll = []
for qind in xrange(0,len(img_names)):
    held = fit_quads[qind]
    if len(betas)!=1:
        mean_beta = np.mean([betas[bind] for bind in xrange(0,len(betas)) if bind != qind],axis=0)
    else:
        mean_beta = betas[0]
    logl = inhomoLogLike(held.pat,mean_beta,temp_folder_name,img_names[qind])
    cv_ll.append(logl)

avg_ll = np.mean(cv_ll)
var_ll = np.var(cv_ll)
with open('likelihood.txt','a') as f:
    f.write('average: '+str(avg_ll)+'\n')
    f.write('variance: '+str(var_ll)+'\n')
    f.write('\n')