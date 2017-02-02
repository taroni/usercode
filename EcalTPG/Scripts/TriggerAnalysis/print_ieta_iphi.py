import ROOT

file0=ROOT.TFile("ECALTPGtree_246926.root")
tree = file0.Get("EcalTPGAnalysis")

entries= tree.GetEntries()
eta4 = [[]]
phi4 = [[]]
evtnb=0
for evt in tree:
    evtnb=evtnb+1
    if (evtnb/100000.) ==  int(evtnb/100000.): print "processing event", evtnb
    for n, flag in enumerate(evt.ttFlag):
        if (evt.ttFlag[n]!=4 ) : continue
        ##print flag, evt.ieta[n], evt.iphi[n]
        #print eta4
        mylist=[]
        matched = False
        matchedIndex = -999
        if eta4[0]!=[]:
            for neta,eta in enumerate(eta4):
                if eta[0]==evt.ieta[n]:
                    matched = True
                    matchedIndex=neta
                    ##print eta[0], evt.ieta[n]
            if matched and matchedIndex!=-999:
                eta4[matchedIndex][1]= 1+eta4[matchedIndex][1]
            else:
                eta4.append( [evt.ieta[n], 1 ] )
        else:
            eta4=[ [evt.ieta[n], 1 ] ]
        if phi4[0]!=[]:
            for nphi,phi in enumerate(phi4):
                if phi[0]==evt.iphi[n]:
                    matched = True
                    matchedIndex=neta
                    ##print phi[0], evt.iphi[n]
            if matched and matchedIndex!=-999:
                phi4[matchedIndex][1]= 1+phi4[matchedIndex][1]
            else:
                phi4.append( [evt.iphi[n], 1 ] )
        else:
            phi4=[ [evt.iphi[n], 1 ] ]
       


myetalist =[eta[0] for eta in eta4]
myphilist= [phi[0] for phi in phi4]
print "phi", myphilist
print "eta", myetalist
