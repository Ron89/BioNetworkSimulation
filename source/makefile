ifeq ($(OUTPUT),0)
	$(OUTPUT) = TEST1
endif

main: 	main.cpp basicDef.h BCNetwork.h gillespie.h  modelLoader.h ODECommon.h ODEOperation.h rungeKutta.h
	g++ -Wall -g main.cpp -o ~/bin/$(OUTPUT) 
