import ROOT
import sys
from array import array
from DataFormats.FWLite import Events, Handle, Runs
import csv

evtList=[]
recHithandle = Handle ("edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >")
recHitLabel = ("ecalRecHit")
recHitProductLabel = ("EcalRecHitsEB")

#events0 = Events("136.85_RunEGamma2018A+RunEGamma2018A+HLTDR2_2018+RECODR2_2018reHLT_skimEGamma_Offline_L1TEgDQM+HARVEST2018_L1TEgDQM_OFF/step3.root")
#events1 = Events("136.85_RunEGamma2018A+RunEGamma2018A+HLTDR2_2018+RECODR2_2018reHLT_skimEGamma_Offline_L1TEgDQM+HARVEST2018_L1TEgDQM/step3.root")

events0 = Events("tmpOFF/step3.root")
events1 = Events("136.85_RunEGamma2018A+RunEGamma2018A+HLTDR2_2018+RECODR2_2018reHLT_skimEGamma_Offline_L1TEgDQM+HARVEST2018_L1TEgDQM_ON_UncalRH/step3_withUncalibRecHit.root")


count=0
oldRun=-99
for event0 in events0:
    ##print event.eventAuxiliary().event(), event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock()
    #Need to pick up the same events in the other file

    event0.getByLabel (recHitLabel, recHitProductLabel, recHithandle)
    recHits0=recHithandle.product()
    if oldRun!=event0.eventAuxiliary().run():
        oldRun=event0.eventAuxiliary().run()
        print 'processing run:', oldRun

    #if count>100: break
    
    events1.to(count)
    #print event0.eventAuxiliary().event(), event0.eventAuxiliary().run(), event0.eventAuxiliary().luminosityBlock()##, event1.eventAuxiliary().event(), event1.eventAuxiliary().run(), event1.eventAuxiliary().luminosityBlock()
    for event1 in events1:
        #print 'event1', event1.eventAuxiliary().event(), event1.eventAuxiliary().run(), event1.eventAuxiliary().luminosityBlock()


        if event1.eventAuxiliary().event()!=event0.eventAuxiliary().event() and event1.eventAuxiliary().run()!=event0.eventAuxiliary().run(): break
        event1.getByLabel (recHitLabel, recHitProductLabel, recHithandle)
        recHits1=recHithandle.product()
        checkRawID=-99
        for n, rechit0 in enumerate(recHits0):
            rechit1=recHits1[n]
            #print rechit0.detid().rawId(), rechit0.chi2(), rechit1.detid().rawId(), rechit1.chi2()
            
            if rechit0.chi2()!=rechit1.chi2():
                evtList.append([event0.eventAuxiliary().run(), event0.eventAuxiliary().luminosityBlock(), event0.eventAuxiliary().event(), rechit0.detid().rawId()])
                print "ATTENTION: chi2 is different (%f, %f) in  run=%d, lumi=%d, evt=%d, detId (%d, %d)" %(rechit0.chi2(), rechit1.chi2(), event0.eventAuxiliary().run(), event0.eventAuxiliary().luminosityBlock(), event0.eventAuxiliary().event(), rechit0.detid().rawId(), rechit1.detid().rawId())
                checkRawID=rechit1.detid().rawId()
                for i in range(0, 19):
                    print i, rechit0.checkFlag(i), rechit1.checkFlag(i)
            if rechit1.detid().rawId()==checkRawID+1 : print rechit0.detid().rawId(), rechit0.chi2(), rechit1.detid().rawId(), rechit1.chi2()
        break
        


    #for recHit in recHits:
    #    print count, recHit.chi2()
    count+=1
with open('evt.csv', 'w') as writeFile:
    writer = csv.writer(writeFile)
    writer.writerows(evtList)
writeFile.close()
