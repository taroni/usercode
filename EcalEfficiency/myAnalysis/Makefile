JOBID         = barrelThreshold90-95
CC            = g++ -g
LIBS          = `root-config --cflags --libs` -L $(ROOFITSYS)/lib -lRooFit -lRooFitCore -I$(ROOFITSYS)/include 

SRCS          = FuncCB.cc fitEfficiency.cc 
SELECTSRC     = selectPairs.cc
MAKESRC       = makePairs.cc

OBJS          = FuncCB.o fitEfficiency.o 

FIT           = fitEfficiency
SELECT        = selectPairs
MAKEPAIRS     = makePairs


all: $(FIT) 

$(FIT):   $(SRCS) $(SELECT) 
	echo "Linking $(FIT) ..."
	g++  -o $(FIT) $(SRCS) $(LIBS)

$(SELECT): $(SELECTSRC) $(MAKEPAIRS) 
	g++  -o $(SELECT) $(SELECTSRC) $(LIBS)
	test -d $(JOBID) || make $(MAKEPAIRS)  make runPairs
	test -d $(JOBID)/selectPairsDir || mkdir $(JOBID)/selectPairsDir

$(MAKEPAIRS): $(MAKESRC) 
	g++  -o $(MAKEPAIRS) $(MAKESRC) $(LIBS)
	test -d $(JOBID) || mkdir $(JOBID)
	test -d $(JOBID)/makePairsDir || mkdir $(JOBID)/makePairsDir


runPairs: $(MAKEPAIRS)
	./$(MAKEPAIRS)  $(JOBID)
runSelect: $(SELECT)
	./$(SELECT) $(JOBID)
runFit: $(FIT)
	./$(FIT) $(JOBID)

runAll: $(FIT)
	./$(MAKEPAIRS)  $(JOBID)
	./$(SELECT) $(JOBID)
	./$(FIT) $(JOBID)


cleanSelect:
	rm -f $(SELECT)

cleanPairs: 
	rm -f $(MAKEPAIRS)
cleanFit:
	rm -f $(FIT)
clean:
	rm -f $(OBJS) core $(FIT) $(SELECT) $(MAKEPAIRS)

