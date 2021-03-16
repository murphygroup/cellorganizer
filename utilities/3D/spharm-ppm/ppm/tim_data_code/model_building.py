from model_building_packages import *

final_output_path=sys.argv[1]
intermediate_output_path=sys.argv[2]
dummy_num=int(sys.argv[3])
rand_num=int(sys.argv[4])
cv_mode=sys.argv[5]
try:
    fold=int(sys.argv[6])
except ValueError:
    if sys.argv[6]=='leave_one_out':
        fold=0
    else:
        print ('Wrong fold parameter')
cv_round=int(sys.argv[7])
mc=sys.argv[8]

ROIs_CenterFilenames=get_files_by_regex(os.path.join(intermediate_output_path,'others'),'centers__.*.txt')
quads = []
for ROI_CenterFilename in ROIs_CenterFilenames:
    pts=np.loadtxt(os.path.join(intermediate_output_path,'others',ROI_CenterFilename), dtype=float)
    img_ROI_name=get_substring_between('centers__','.txt',ROI_CenterFilename)
    maxx,maxy,maxz=find_xmax_ymax_zmax(intermediate_output_path,img_ROI_name)
    # print maxx,maxy,maxz
    minx=1
    miny=1
    minz=1
    pts_scale=np.column_stack(((pts[:,0]-minx)/(maxx-minx),(pts[:,1]-miny)/(maxy-miny),(pts[:,2]-minz)/(maxz-minz)))
    cur_pat=PointPattern(pts_scale,dimension=3)
    print (ROI_CenterFilename+(' PointPattern complete'))
    quads.append(QuadScheme(cur_pat,dummy_num))
    print (ROI_CenterFilename+(' QuadScheme complete'))

img_roi_names=[get_substring_between('centers__','.txt',ROI_CenterFilename) for ROI_CenterFilename in ROIs_CenterFilenames]
full_model,betas,failed=FitModel_all(quads,intermediate_output_path,img_roi_names)
with open(os.path.join(final_output_path,'full_model.txt'),'a') as f:
    for b in full_model:
        f.write(str(b)+' ')
    f.write('\n')

#save full model as mat
factor_channel_num=get_factor_channel_num()
fcs=''
for num in factor_channel_num:
    fcs=fcs+str(num)+'-'
fcs=fcs[0:-1]
sio.savemat(os.path.join(final_output_path,'model__mc-'+mc+'fc-'+fcs+'.mat'),{'coefficients': full_model})

# fit_quads = [quad for qind,quad in enumerate(quads) if qind not in failed]


img_names=set2list([get_substring_between('centers__','-[0-9]+.txt',ROI_CenterFilename) for ROI_CenterFilename in ROIs_CenterFilenames])
if cv_mode=='rd_img':
    element_names=img_names
if cv_mode=='rd_roi':
    element_names=img_roi_names

if fold==0:
    fold==len(element_names)

if len(betas)!=1:
    element_ids=[rd for rd in range(0,len(element_names))]
    group_size=len(element_names)//fold
    cv_ll=[]
    #loop cv_round
    for cv_r in range(0,cv_round):
        random.shuffle(element_ids)
        #get cv_assignment
        cv_assignment=[]
        for fold_ind in range(0,fold):
            cv_assignment.append([fold_ind*group_size+group_ind for group_ind in range(0,group_size)])
        if cv_mode=='rd_img':
            cv_assignment=get_cv_assignment_real(cv_assignment,img_names,img_roi_names)
        #loop fold
        if len(cv_assignment)!=1:
            for fold_ind in range(0,fold):
                # if cv_mode=='rd_img':
                #     test_img_roi_names=cv_assignment[fold_ind]
                #     ["foo", "bar", "baz"].index(test_img_roi_name) for test_img_roi_name in test_img_roi_names
                # if cv_mode=='rd_roi':
                test_set=cv_assignment[fold_ind]
                test_img_roi_names=[img_roi_names[ind] for ind in test_set]
                train_set=[cv_assignment[cv_a_ind] for cv_a_ind in range(0,len(cv_assignment)) if cv_a_ind!=fold_ind]
                print (test_img_roi_names)
                print (train_set)

                mean_betas=[]
                for group in train_set:
                    mean_betas.append(np.mean([betas[element] for element in group],axis=0))
                mean_beta=np.mean(mean_betas,axis=0)

                logl=inhomoLogLike(test_img_roi_names,mean_beta,intermediate_output_path,rand_num)
                cv_ll.append(logl)
        else:
            test_img_roi_names=[img_roi_names[ind] for ind in cv_assignment[0]]
            mean_beta=np.mean([betas[element] for element in range(0,len(test_img_roi_names))],axis=0)
            logl=inhomoLogLike(test_img_roi_names,mean_beta,intermediate_output_path,rand_num) ###test=train
            cv_ll.append(logl)

    # cv_ll = []
    # for qind in range(0,len(img_names)):
    #     held = fit_quads[qind]
    #     if len(betas)!=1:
    #         mean_beta = np.mean([betas[bind] for bind in range(0,len(betas)) if bind != qind],axis=0)
    #     else:
    #         mean_beta = betas[0]
    #     logl = inhomoLogLike(held.pat,mean_beta,temp_folder_name,img_names[qind])
    #     cv_ll.append(logl)

    avg_ll = np.mean(cv_ll)
    var_ll = np.var(cv_ll)
    with open(os.path.join(final_output_path,'likelihood.txt'),'a') as f:
        f.write('average: '+str(avg_ll)+'\n')
        f.write('variance: '+str(var_ll)+'\n')
        f.write('\n')
else:
    # mean_beta=betas[0]
    # cv_ll=[inhomoLogLike(test_img_roi_names,mean_beta,intermediate_output_path,rand_num)] ###test=train
    print ('Only one image and one ROI, cannot get likelihood. But the model was still saved.')
