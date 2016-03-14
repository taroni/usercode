import ROOT
import math
import numpy as np
hfile = ROOT.TFile("TheoUnc_HiggsToMuTau_125_VBF_SM.root", "READ")
hdir = hfile.Get("demo") 
hlist = hdir.GetListOfKeys()

mupt_lep = [50,45,25]
mupt_had = [45,35,30]
ept = [15,15,15]
taupt = [35,40,40]

hWeights=[]
weights = []

histos=[]
histosAcc=[]
refvaluennpdf=[]
maxnnpdf=[-10., -10.,-10.]
minnnpdf=[100., 100., 100.]

alphannpdf=[]


for j in range(0, 3):
    histo=ROOT.TH1F('histoNNPDF_'+str(j), 'histoNNPDF_'+str(j), 3000, 0.,0.0003)
    histoAcc=ROOT.TH1F('histoAccNNPDF_'+str(j), 'histoAccNNPDF_'+str(j), 3000, 0.,0.03)
    histos.append(histo)
    histosAcc.append(histoAcc)
for hname in hlist:
    if 'Weights' in hname.GetName():
        # print hname
        hWeights.append(hname.GetName())
        weights.append((hname.GetName(), hdir.Get(hname.GetName()).Integral()))

        
        
for j in range(0, 3):
    had_pdf=[]
    lep_pdf=[]
    had_scale=[]
    lep_scale=[]
    h_had_pdf_2k =[]
    h_had_pdf_3k =[]
    h_had_pdf_4k =[]
    h_hadscale =[]
    h_lep_pdf_2k =[]
    h_lep_pdf_3k =[]
    h_lep_pdf_4k =[]
    h_lepscale =[]
    for hname in hlist:
        
        if hdir.Get(hname.GetName()).IsA() != ROOT.TH2F.Class() : continue
        if str(j)+'j' not in hname.GetName(): continue
        
        if 'had' in hname.GetName():
            if 'j40' in hname.GetName():
                h_had_pdf_4k.append(hname.GetName())
            if 'j30' in hname.GetName():
                h_had_pdf_3k.append(hname.GetName())
            if 'j2' in hname.GetName():
                h_had_pdf_2k.append(hname.GetName())
            if 'j10' in hname.GetName():
                h_hadscale.append(hname.GetName())
        if 'lep' in hname.GetName():
            if 'j40' in hname.GetName():
                h_lep_pdf_4k.append(hname.GetName())
            if 'j30' in hname.GetName():
                h_lep_pdf_3k.append(hname.GetName())
            if 'j2' in hname.GetName():
                h_lep_pdf_2k.append(hname.GetName())
            if 'j10' in hname.GetName():
                h_lepscale.append(hname.GetName())


    
    ref_had_scale =0.
    ref_lep_scale =0.
    ##had scale syst
    int_had_scale=[]
    int_lep_scale=[]
    totalIntegral = []
    totalEntries = []

    
    ref_had_scale = 0
    ref_lep_scale=0
    
    for n,name in enumerate(h_hadscale):
        
        if '1001' in name : 
            h_ref_had_scale = hdir.Get(name)
            ref_had_scale = h_ref_had_scale.Integral(taupt[j]+1, 100, mupt_had[j]+1, 100)/weights[0][1]
            h_ref_lep_scale = hdir.Get(h_lepscale[n])
            ref_lep_scale = h_ref_lep_scale.Integral(ept[j]+1, 100, mupt_lep[j]+1, 100)/weights[0][1]
            totalIntegral.append(h_ref_had_scale.Integral())
            totalEntries.append(h_ref_had_scale.GetEntries())
        else:
            h_ref_had_scale = hdir.Get(name)
            h_ref_lep_scale = hdir.Get(h_lepscale[n])
            int_had_scale.append( h_ref_had_scale.Integral(taupt[j]+1, 100, mupt_had[j]+1, 100)/weights[n][1])
            int_lep_scale.append( h_ref_lep_scale.Integral(ept[j]+1, 100, mupt_lep[j]+1, 100)/weights[n][1])
            totalIntegral.append(h_ref_had_scale.Integral())
            totalEntries.append(h_ref_had_scale.GetEntries())
            
    for n,value in enumerate(int_had_scale):
        int_had_scale[n]=value+int_lep_scale[n]
    ref_had_scale+=ref_lep_scale
    
    unc= (max(int_had_scale)-min(int_had_scale))/2/ref_had_scale
    
    #print (min(int_had_scale)-ref_had_scale)/ref_had_scale,(max(int_had_scale)-ref_had_scale)/ref_had_scale
    print 'scale %s jet: %.1f or %.1f/%.1f + Add yellow report reccomendation' %(str(j),100*unc, 100*(min(int_had_scale)-ref_had_scale)/ref_had_scale,100*(max(int_had_scale)-ref_had_scale)/ref_had_scale)
    
    
    ref_had_scale =0.
    ref_lep_scale =0.
    int_had_scale=[]
    int_lep_scale=[]
    int_pdf=[]

    lower=[]
    higher=[]
    
    h_ref_had_scale = hdir.Get(h_hadscale[0])
    ref_had_scale = h_ref_had_scale.Integral(taupt[j]+1, 100, mupt_had[j]+1, 100)/weights[0][1]
    h_ref_lep_scale = hdir.Get(h_lepscale[0])
    ref_lep_scale = h_ref_lep_scale.Integral(ept[j]+1, 100, mupt_lep[j]+1, 100)/weights[0][1]
    nnpdf=0
    for k,w in enumerate(weights):
        #print k, w
        if '2001' in w[0]:
            nnpdf=k
    for n,name in enumerate(h_had_pdf_2k):
        h_ref_had_scale = hdir.Get(name)
        h_ref_lep_scale = hdir.Get(h_lep_pdf_2k[n])
        int_had_scale.append( h_ref_had_scale.Integral(taupt[j]+1, 100, mupt_had[j]+1, 100)/weights[nnpdf+n][1])
        int_lep_scale.append( h_ref_lep_scale.Integral(ept[j]+1, 100, mupt_lep[j]+1, 100)/weights[nnpdf+n][1])

    for n,value in enumerate(int_had_scale):
        int_had_scale[n]=value+int_lep_scale[n]
        
    ref_had_scale+=ref_lep_scale
    

    
    a_nnpdf=ref_had_scale
    refvaluennpdf.append(a_nnpdf)    

    average=np.mean(int_had_scale)
    #print ref_had_scale, int_had_scale, average   
    
    for n in range(0, len(int_had_scale)-2):
        histosAcc[j].Fill(int_had_scale[n])
        histos[j].Fill(int_had_scale[n]-ref_had_scale)
        #print int_had_scale[n]-ref_had_scale
        if (int_had_scale[n]-ref_had_scale)>maxnnpdf[j]:
            maxnnpdf[j]= int_had_scale[n]-ref_had_scale
        if(int_had_scale[n]-ref_had_scale)<minnnpdf[j]:
            minnnpdf[j]= int_had_scale[n]-ref_had_scale


        
 

        if (int_had_scale[n]-average)> 0 :
            higher.append(int_had_scale[n])
        else:
            lower.append(int_had_scale[n])


    if len(higher)!=0 : maxhigh=max(higher)
    if len(lower)!=0 :  minlow= min(lower)

    for n in range (1, 10):
        diff= abs(minlow-average)
        thr = minlow+n*diff/10.
        if len(lower)==0 : break
        count=len([x for x in lower if x>thr])
        if count<0.68*(len(lower)+1) :
            print 'unc %.0f : %.3f' %(j, 100*0.5*(minlow+((n-1)*diff/10. - n*diff/10.)-average)/average)
            break
        
    for n in range (1, 10):
        diff= abs(maxhigh-average)
        thr = maxhigh-n*diff/10.
        if len(higher)==0 : break
        count=len([x for x in higher if x<thr ])
        if count<0.68*(len(higher)+1) :
            print 'unc %.0f : %.3f' %(j, 100*0.5*(maxhigh-(n*diff/10. - (n-1)*diff/10.)-average)/average)
            break
    
    listsize=len(int_had_scale)
    print int_had_scale[listsize-2],int_had_scale[listsize-1], int_had_scale[listsize-2]-int_had_scale[listsize-1]
    alphannpdf.append(1.333*(int_had_scale[listsize-2]-int_had_scale[listsize-1])/2)

 


    
        
outfile = ROOT.TFile("nnpdf_"+hfile.GetName(), "RECREATE")
outfile.cd()
for j,h in enumerate(histos):
    rms=h.GetRMS()
    print 'nnpdf unc for %.0f jet: %2.1f or %2.1f/%2.1f + alpha %2.2f' %(j, 100*rms/refvaluennpdf[j], 100*minnnpdf[j]/refvaluennpdf[j],  100*maxnnpdf[j]/refvaluennpdf[j], 100*alphannpdf[j]/refvaluennpdf[j])
    h.Write()
    histosAcc[j].Write()


    
outfile.Close()


