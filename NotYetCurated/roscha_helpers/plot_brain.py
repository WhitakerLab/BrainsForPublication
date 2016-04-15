#!/usr/bin/env python2.7
def plotZslices_alloption(niftipath,mnipath='',ortho='z',cut_coords='',Nraw=1,smoothing=0,LR=False,outdir='',colorpos='r',colorneg='b',Zannotate=False,thresholdpos='def',Zannotates='def',thresholdneg=False,alphamap=1,alphabrain=1):
    "niftipath: path to the nifti file, can be a 3D - if activation map, specify thresholds,"
    "mnipath : path to the mni T1 brain
    "cut_coords can be a int as the number of zslices to display of a list of slices number (in MNI) (even list of one to get one specific slice)"
    "Nraw: the number of raw"
    "smoothing: number of voxel to smooth; LR:annotate left and right"
    "outdir:path to save the file"
    "color:list of color for each volume, or only one color, neg or pos if corresponding threshold to display"
    "Zannotate : Number=annotate z number, False=not annotate, Brain=on a X slice, with lign, or Both"
    "thresholdpos: specify threshold to cut and see above (can be a list for activation map: layer effect) or False will not be displayed or 'def' as 0.5 on normalized file"
    "thresholdneg: specify threshold to cut and see bellow (can be a list for activation map: layer effect) or False will not be displayed "
    import matplotlib.pyplot as plt
    import numpy as np
    import nilearn.plotting
    import nilearn.image   
    from nilearn.plotting import plot_roi, plot_stat_map
    from nilearn.plotting.find_cuts import find_cut_slices
    from nilearn.image.image import mean_img
    import nibabel
    import seaborn as sns
    initialcol=sns.light_palette((0,0,0), as_cmap=True)#'Greys'
    data=nibabel.load(niftipath)
    datasize=data.get_shape()
    lineW=1./(Nraw+int((Zannotate=='Brain' or Zannotate=='Both')))
    if mnipath=='':
        mnipath='/home/mrstats/maamen/Software/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz' ##this only works for the donders institute (Nijmegen, The Neterlands)
    if type(cut_coords)==int or cut_coords=='':
        if cut_coords=='':
            cut_coords=6
        #find best cut
        if len(datasize)==4:
            #for 4D nifti
            cut_coords=find_cut_slices(mean_img(nibabel.nifti1.Nifti1Image(np.sign(np.abs(data.get_data())),data.get_affine())), n_cuts=cut_coords)
        else:
            #for 3D nifti
            cut_coords=find_cut_slices(data, n_cuts=cut_coords)
    
    #split in N raw
    if cut_coords!=(0,0,0):
        cut_coords=np.array(cut_coords)
        cc=cut_coords
        cut_coords=[cut_coords[i*len(cut_coords)/np.float(Nraw):(i+1)*len(cut_coords)/np.float(Nraw)] for i in range(Nraw)]
    else:
        cut_coords=[cut_coords]
    #define color as a vector (length :the number of volume):
    #if not enought color are proveded, the last of them is repeated
    #color are defined independantly for negative value display and positive value display
    if type(colorneg)==str:
        colorneg=[colorneg]
    if type(colorpos)==str:
        colorpos=[colorpos]
    if len(datasize)==4 and len(colorpos)!=datasize[3]:
        provcol=colorpos[len(colorpos)-1]
        colorpos=np.concatenate([colorpos,[provcol for i in range(datasize[3]-len(colorpos))]])
    if len(datasize)==4 and len(colorneg)!=datasize[3]:
        provcol=colorneg[len(colorneg)-1]
        colorneg=np.concatenate([colorneg,[provcol for i in range(datasize[3]-len(colorneg))]])

    #adjust threshold by normalizing image in the default version and taking 0.5
    if thresholdpos=='def':
        data=nibabel.nifti1.Nifti1Image(data.get_data()/np.float(np.max(data.get_data())),data.get_affine())
        thresholdpos=[0.5]
    #organize thresholds, more than 1 threshold to make a layer effect,
    #positive and negative values display are treated independantly:
    if type(thresholdpos)!=np.bool: thresholdpos=[i for i in np.sort(thresholdpos)]
    if type(thresholdneg)!=np.bool: thresholdneg=[i for i in -np.sort(-1*np.array(thresholdneg))]
    
    #load data to create a white backgroung
    func=mean_img(nibabel.nifti1.Nifti1Image(np.sign(np.abs(data.get_data())),data.get_affine()))
    
    ####################subplot
    for i in range(Nraw):
        
        ax=plt.subplot(Nraw+int((Zannotate=='Brain' or Zannotate=='Both')),1,i+1)
        
        #plot the white backgroung as a zeros value brain (without it, the view focus aroung the first area plotted)
        brain=nilearn.plotting.plot_roi(nibabel.nifti1.Nifti1Image(np.zeros(func.get_shape()),data.get_affine()), nibabel.nifti1.Nifti1Image(np.zeros(func.get_shape()),data.get_affine()),colorbar=False,cut_coords=cut_coords[i],display_mode=ortho,alpha=1,draw_cross=False,cmap=initialcol,black_bg=False,axes=ax,annotate=False)
        
        ###############plot the volumes for Z brain slices
        if len(datasize)==3:
            iter_imgs=[data]
        else:
            iter_imgs=nilearn.image.iter_img(niftipath)
        
        j=0
        for img in iter_imgs:
            
            ##plot the positive values
            if thresholdpos!=False:
                colorprovpos=sns.light_palette(colorpos[j],len(thresholdpos),reverse=True)[::-1]
            
                
                img2=nilearn.image.smooth_img(img,smoothing)
                ##plot the different threshold (layer effect) for the positive values
                for kn,k in enumerate(thresholdpos):
                    brain.add_contours(img2,filled=True,levels=[k],cmap=None,colors=[[colorprovpos[kn][0],colorprovpos[kn][1],colorprovpos[kn][2]]],linewidths=lineW,alpha=k/np.max(thresholdpos))#alphamap)                    
            ##plot the negative values   
            if thresholdneg!=False:
                colorprovneg=sns.light_palette(colorneg[j],len(thresholdneg),reverse=True)[::-1]
                #switch negative to positive
                img2=nibabel.nifti1.Nifti1Image(-1*nilearn.image.smooth_img(img,smoothing).get_data(),data.get_affine())
                ##plot the negatives values for each negative threshold
                for kn,k in enumerate(thresholdneg):
                    brain.add_contours(img2,filled=True,levels=[k],cmap=None,colors=[[colorprovneg[kn][0],colorprovneg[kn][1],colorprovneg[kn][2]]],linewidths=lineW,alpha=k/np.max(thresholdpos))#alphamap)                    
                    
            j+=1
        ##plot the brain contour for Z brain slices
        #externe
        brain.add_contours(nilearn.image.smooth_img(mnipath,5),alpha=1*alphabrain, levels=[95],linewidths=lineW, cmap=sns.dark_palette('w', as_cmap=True),)       
        #interne (a little transparent)
        brain.add_contours(nilearn.image.smooth_img(mnipath,0.5),alpha=0.8*alphabrain, levels=[5000],linewidths=lineW)
        #add annotation if reauested
        if Zannotate=='Both' or Zannotate=='Number' :
            brain.annotate(left_right=LR,size=int(12*lineW))
        
        print 'raw '+str(i)+' ready'

        
    ########################## plot the X brain (same process but on X)
    if Zannotate=='Brain' or Zannotate=='Both':
        print 'doing annotate X slice'
        ax=plt.subplot(Nraw+1,1,Nraw+1)
        if len(datasize)==4:
            cut_coords=find_cut_slices(mean_img(nibabel.nifti1.Nifti1Image(np.sign(np.abs(data.get_data())),data.get_affine())), n_cuts=1,direction='x')
        else:
            cut_coords=find_cut_slices(data, n_cuts=1,direction='x')
        
        #plot the white background Xbrain
        brain=nilearn.plotting.plot_roi(nibabel.nifti1.Nifti1Image(np.zeros(func.get_shape()),data.get_affine()), nibabel.nifti1.Nifti1Image(np.zeros(func.get_shape()),data.get_affine()),colorbar=False,cut_coords=cut_coords,display_mode='x',alpha=1,draw_cross=False,cmap=initialcol,black_bg=False,axes=ax,annotate=False)
        
        if Zannotate=='Both' or Zannotate=='Number' :
            brain.annotate(left_right=LR,size=int(12*lineW))
            
        if len(datasize)==3:
            iter_imgs=[data]
        else:
            iter_imgs=nilearn.image.iter_img(niftipath)
        #plot the volumes
        j=0
        for img in iter_imgs:
            if thresholdpos!=False:
                colorprovpos=sns.light_palette(colorpos[j],len(thresholdpos),reverse=True)[::-1]
                img2=nilearn.image.smooth_img(img,smoothing)
                for kn,k in enumerate(thresholdpos):
                    brain.add_contours(img2,filled=True,levels=[k],cmap=None,colors=[[colorprovpos[kn][0],colorprovpos[kn][1],colorprovpos[kn][2]]],linewidths=lineW,alpha=k/np.max(thresholdpos))#alphamap)                    
                    
            if thresholdneg!=False:
                colorprovneg=sns.light_palette(colorneg[j],len(thresholdneg),reverse=True)[::-1]
                img2=nibabel.nifti1.Nifti1Image(-1*nilearn.image.smooth_img(img,smoothing).get_data(),data.get_affine())
                for kn,k in enumerate(thresholdneg):
                    brain.add_contours(img2,filled=True,levels=[k],cmap=None,colors=[[colorprovneg[kn][0],colorprovneg[kn][1],colorprovneg[kn][2]]],linewidths=lineW,alpha=k/np.max(thresholdpos))#alphamap)                    
                
            j+=1
        brain.add_contours(nilearn.image.smooth_img(mnipath,5),alpha=1*alphabrain, levels=[95],linewidths=lineW, cmap=sns.dark_palette('w', as_cmap=True),)       
        brain.add_contours(nilearn.image.smooth_img(mnipath,0.5),alpha=0.8*alphabrain, levels=[5000],linewidths=lineW)
        ##plot the line indicating the cut
        for i in cc:
            ax.plot([-100, 100], [i, i], 'k-',lw=lineW)#/(85.+73.)
        ax.axis((-300.0, 300.0, -80.0, 110.0))   
    #save
    if outdir!='':
        plt.savefig(outdir,dpi=300)

