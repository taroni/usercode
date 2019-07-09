import ROOT
import sys
from array import array
from DataFormats.FWLite import Events, Handle, Runs
import csv

evtList=[]
with open('evt.csv', 'r') as readFile:
    reader = csv.reader(readFile)
    evtList.extend(list(reader))

#print evtList
newEvtList=[]

recHithandle = Handle ("edm::SortedCollection<EcalUncalibratedRecHit,edm::StrictWeakOrdering<EcalUncalibratedRecHit> >")
recHitLabel = ("ecalMultiFitUncalibRecHit")
recHitProductLabel = ("EcalUncalibRecHitsEB")

events0 = Events("tmpOFF/step3.root")
events1 = Events("136.85_RunEGamma2018A+RunEGamma2018A+HLTDR2_2018+RECODR2_2018reHLT_skimEGamma_Offline_L1TEgDQM+HARVEST2018_L1TEgDQM_ON_UncalRH/step3_withUncalibRecHit.root")

#chStatus= Handle("edm::ESHandle<EcalChannelStatus>")

#es.get<EcalChannelStatusRcd>().get(chStatus);
#EcalChannelStatusMap::const_iterator chit = chStatus->find(detid);
#EcalChannelStatusCode::Code dbstatus = chit->getStatusCode();
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
            ##print rechit0.id().rawId(), rechit0.chi2(), rechit1.id().rawId(), rechit1.chi2()
            if [str(event0.eventAuxiliary().run()), str(event0.eventAuxiliary().luminosityBlock()), str(event0.eventAuxiliary().event()), str(rechit0.id().rawId())] in evtList:
                newEvtList.append([str(event0.eventAuxiliary().run()), str(event0.eventAuxiliary().luminosityBlock()), str(event0.eventAuxiliary().event()), str(rechit0.id().rawId())])
                print "CHECKING CHI2  (%f, %f), amplitudes (%f, %f) in  run=%d, lumi=%d, evt=%d, detId (%d, %d)" %(rechit0.chi2(), rechit1.chi2(),rechit0.amplitude(), rechit1.amplitude(), event0.eventAuxiliary().run(), event0.eventAuxiliary().luminosityBlock(), event0.eventAuxiliary().event(), rechit0.id().rawId(), rechit1.id().rawId())

            if rechit0.chi2()!=rechit1.chi2():
                if [str(event0.eventAuxiliary().run()), str(event0.eventAuxiliary().luminosityBlock()), str(event0.eventAuxiliary().event()), str(rechit0.id().rawId())] in evtList:
                    print "ATTENTION: chi2 is different (%f, %f), amplitudes (%f, %f) in  run=%d, lumi=%d, evt=%d, detId (%d, %d)" %(rechit0.chi2(), rechit1.chi2(), rechit0.amplitude(), rechit1.amplitude(),event0.eventAuxiliary().run(), event0.eventAuxiliary().luminosityBlock(), event0.eventAuxiliary().event(), rechit0.id().rawId(), rechit1.id().rawId())
                    
                    ##        #for i in range(0, 19):
                    ##        #    print i, rechit0.checkFlag(i), rechit1.checkFlag(i)
                    ##
                #else:
                #    print "ATTENTION: chi2 is different (%f, %f) but not in the evt list," %(rechit0.chi2(), rechit1.chi2()), [str(event0.eventAuxiliary().run()), str(event0.eventAuxiliary().luminosityBlock()), str(event0.eventAuxiliary().event()), str(rechit0.id().rawId())] 
            ####if rechit1.id().rawId()==checkRawID+1 : print rechit0.id().rawId(), rechit0.chi2(), rechit1.id().rawId(), rechit1.chi2()
        break
        


    count+=1

#print newEvtList

#diffList = [x for x in evtList if x not in newEvtList]
#print diffList
