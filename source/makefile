ifeq ($(OUTPUT),0)
	$(OUTPUT) = TEST1
endif



main: main.cpp histoGen.h modelLoader.h gillespie.h basicDef.h timeDerivative.h
	g++ -Wall -g main.cpp -o ~/bin/$(OUTPUT) 
