ifeq ($(OUTPUT),0)
	$(OUTPUT) = TEST1
endif

main: 	main.cpp basicDef.h BCNetwork.h gillespie.h  modelLoader.h ODECommon.h ODEOperation.h rungeKutta.h marginalDistriExtractor.h
	g++44 -Wall -g -std=c++0x main.cpp -o ~/bin/$(OUTPUT) 
