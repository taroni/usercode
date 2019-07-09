import sys
import re
import csv

logfile=open("136.85_RunEGamma2018A+RunEGamma2018A+HLTDR2_2018+RECODR2_2018reHLT_skimEGamma_Offline_L1TEgDQM+HARVEST2018_L1TEgDQM/step3_RunEGamma2018A+RunEGamma2018A+HLTDR2_2018+RECODR2_2018reHLT_skimEGamma_Offline_L1TEgDQM+HARVEST2018_L1TEgDQM.log", "r")

lines= [line.rstrip() for line in logfile]
interestingLines=[]
newLine=[]
for n,line in enumerate(lines):
    if 'Begin processing' in line:
        interestingLines.append(tuple(newLine))
        del newLine[:]
        newLine.append(line)
    elif 'recovery not possible:' in line:
        newLine.append(line)
    if n==len(lines)-1:
        interestingLines.append(tuple(newLine))

del(lines)

#print interestingLines


interestingLines=[x for x in interestingLines if len(x)>1]

textlist=[]

for line in interestingLines:
    a = re.search("Run", line[0])   
    b = re.search(", Event", line[0])
    c = re.search(", LumiSection", line[0])
    d = re.search(" on stream", line[0])
    
    evt=[line[0][a.end():b.start()].strip(), line[0][c.end():d.start()].strip(), line[0][b.end():c.start()].strip()]
    
    ncrystals= len(line)-1
    #print 'number of crystals', ncrystals
    for nxtal in range(ncrystals):
        nxtal=nxtal+1
        e=re.search(', id=', line[nxtal])
        
        xtalId=line[nxtal][e.end():].strip()
        reason=line[nxtal][:e.start()].strip()
        
        #print 'read:', xtalId, reason
        
        tmplist=evt[:]
        #print evt, tmplist
        tmplist.extend([xtalId, reason])
        #print 'tmplist', tmplist
        textlist.append(tmplist)
        #print 'last added', textlist[-1]



evtList=[]
with open('evt.csv', 'r') as readFile:
    reader = csv.reader(readFile)
    evtList.extend(list(reader))

failedRecov={}
for line in textlist:
    evt=(line[0], line[1], line[2], line[3])
    failedRecov[evt]=line[4]
for line in evtList:
    if tuple(line) in failedRecov.keys():
        print 'unrecovered Event', line, failedRecov[tuple(line)]
    else:
        print 'unknown reason for', line
logfile.close()
