# compose reaction network
# module created for no. of mole. description of network.
# required packages
import copy
import numpy as np
import random
import time
import os
import pickle

class reactionNetwork:
    __rate_size__=4
    __depend_size__=3
    __upd_size__=3
    __index_length__=6
    reactant={}
    reaction={}
    
    def __init__(self):
        random.seed(time.gmtime())
    
    def insert_new_reactant(self,reactant_alias,unexpected=0):
        if unexpected ==1:
            print "reactant",reactant_alias,"not exists, creating instance."
        self.reactant[hash(reactant_alias)]=[reactant_alias,int(raw_input("type in initial condition:"))]
        # (above) change int into float for concentration description
        if unexpected ==1:
            return hash(reactant_alias)
    
    def get_reactant_hash(self,reactant_alias):
        reactant_hash=hash(reactant_alias)
        if reactant_hash in self.reactant:
            return reactant_hash
        else:
            return self.insert_new_reactant(reactant_alias,1)
    
    def list_reactant(self):
        print "format:\n index, [reactant name, initial condition]"
        for i in self.reactant:
            print i, self.reactant[i]
    
    def insert_new_reaction(self,code_alias, rate_alias, depend_alias, upd_no_alias, upd_alias):
    # insert_new_reaction(int, list, list, int, list)
        # generate index for current reaction
        reaction_index=random.randint(10**self.__index_length__,10**(self.__index_length__+1)-1)
        while reaction_index in self.reaction:
            reaction_index=random.randint(10**self.__index_length__,10**(self.__index_length__+1)-1)
        # generate temporary structure to store reaction
        temp_react=[code_alias,np.zeros(self.__rate_size__),np.zeros(self.__depend_size__, dtype=int),upd_no_alias,[]]
        for i in range(self.__upd_size__):
            temp_react[-1].append([0,0])
        # insert reaction information
        if code_alias==0:
            temp_react[1][0] = rate_alias[0]
        elif code_alias==1:
            temp_react[1][0] = rate_alias[0]
            temp_react[2][0] = self.get_reactant_hash(depend_alias[0])
        elif code_alias==2:
            temp_react[1][0] = rate_alias[0]
            temp_react[2][0] = self.get_reactant_hash(depend_alias[0])
            temp_react[2][1] = self.get_reactant_hash(depend_alias[1])
        elif code_alias==3:
            temp_react[1][0] = rate_alias[0]
            temp_react[1][1] = rate_alias[1]
            temp_react[1][2] = rate_alias[2]
            temp_react[2][0] = self.get_reactant_hash(depend_alias[0])
            temp_react[2][1] = self.get_reactant_hash(depend_alias[1])
        elif code_alias==4:
            temp_react[1][0] = rate_alias[0]
            temp_react[1][1] = rate_alias[1]
            temp_react[2][0] = self.get_reactant_hash(depend_alias[0])
            temp_react[2][1] = self.get_reactant_hash(depend_alias[1])
        else:
            #extend cases here
            print "unknown reaction type, fatal",
            return
        if upd_no_alias<=self.__upd_size__:
            for i in range(upd_no_alias):
                temp_react[4][i][0]=self.get_reactant_hash(upd_alias[i][0])
                temp_react[4][i][1]=upd_alias[i][1]
        else:
            print "exceed maximum updated reactants, fatal"
            return
        self.reaction[reaction_index]=copy.deepcopy(temp_react)
        
    def list_reaction(self):
        print "format:\n index: reaction type\t reaction rate \n reaction\t reaction rate constants\n"
        for i in self.reaction:
            # generate RHS of equation
            RHS_reactant={}
            for j in self.reactant:
                RHS_reactant[j]=0
            print i,':',
            if self.reaction[i][0]==0:
                print "order-0\t",
                print "k =", self.reaction[i][1][0]
                print "Nul. ->",
            elif self.reaction[i][0]==1:
                print "order-1\t",
                print "k =", self.reaction[i][1][0]
                print self.reactant[self.reaction[i][2][0]][0],"->",
                RHS_reactant[self.reaction[i][2][0]]+=1
            elif self.reaction[i][0]==2:
                print "order-2\t",
                print "k =", self.reaction[i][1][0]
                print self.reactant[self.reaction[i][2][0]][0],'+',self.reactant[self.reaction[i][2][1]][0],"->",
                RHS_reactant[self.reaction[i][2][0]]+=1
                RHS_reactant[self.reaction[i][2][1]]+=1
            elif self.reaction[i][0]==3:
                print "Hill equation\t",
                print "k =", self.reaction[i][1][0], "\tK =", self.reaction[i][1][1], "\tH =", self.reaction[i][1][2]
                print self.reactant[self.reaction[i][2][0]][0],'+',self.reaction[i][1][2],self.reactant[self.reaction[i][2][1]][0],"<=> [complex] ->",
                RHS_reactant[self.reaction[i][2][0]]+=1
                RHS_reactant[self.reaction[i][2][1]]+=self.reaction[i][1][2]
            elif self.reaction[i][0]==4:
                print "MM kinetics\t",
                print "k =", self.reaction[i][1][0], "\tK =", self.reaction[i][1][1]
                print self.reactant[self.reaction[i][2][0]][0],'+',self.reactant[self.reaction[i][2][1]][0],"<=> [complex] ->",
                RHS_reactant[self.reaction[i][2][0]]+=1
                RHS_reactant[self.reaction[i][2][1]]+=1
            else:
                # extend cases here
                print "undefined case, fatal"
                return
            for j in range(self.reaction[i][3]):
                RHS_reactant[self.reaction[i][4][j][0]]+=self.reaction[i][4][j][1]
            tempMarker=0
            for j in RHS_reactant:
                if RHS_reactant[j]!=0:
                    if tempMarker==0:
                        tempMarker=1
                    else:
                        print '+',
                    if RHS_reactant[j]==1:
                        print self.reactant[j][0],
                    else:
                        print RHS_reactant[j],self.reactant[j][0],
            if tempMarker==0:
                print 'Nul.',
            print "\n"
    def export(self,exportName):
        if os.path.exists("./"+exportName)==0:
            os.mkdir(exportName)
        reactant_export=copy.deepcopy(self.reactant)
        reaction_export=copy.deepcopy(self.reaction)
        index_temp=0
        export_readme=open(exportName+'/readme','w')
        export_readme.write("Exported time: %2d/%2d/%4d.\n\n"%(time.localtime().tm_mday,time.localtime().tm_mon,time.localtime().tm_year))
        export_readme.write("Reactants:\n")
        export_file=open(exportName+'/reactant','w')
        for i in reactant_export:
            reactant_export[i].append(index_temp)
            index_temp+=1
            export_file.write('%d\n'%reactant_export[i][1])
            export_readme.write('%-2d\t%-8s\t%d\n'%(reactant_export[i][-1],reactant_export[i][0],reactant_export[i][1]))
        export_file.close()
        export_readme.write("\nReactions\n")
        export_file=open(exportName+'/reaction','w')
        index_temp=0
        for i in reaction_export:
            for j in range(min(reaction_export[i][0],2)):
                reaction_export[i][2][j]=reactant_export[reaction_export[i][2][j]][-1]
            for j in range(reaction_export[i][3]):
                reaction_export[i][4][j][0]=reactant_export[reaction_export[i][4][j][0]][-1]
            export_readme.write("%-4d\t"%(index_temp))
            index_temp+=1
            export_file.write('%1d\t'%reaction_export[i][0])
            export_readme.write('Type: %1d\t'%reaction_export[i][0])
            for j in range(self.__rate_size__):
                export_file.write('%.4e\t'%reaction_export[i][1][j])
#                export_readme.write('%.4e\t'%reaction_export[i][1][j])
            export_readme.write('Dependant: ')
            for j in range(self.__depend_size__):
                export_file.write('%-2d\t'%reaction_export[i][2][j])
                if self.reaction[i][2][j]!=0:
                    export_readme.write('%-8s\t'%(self.reactant[self.reaction[i][2][j]][0]))
            export_file.write('%-2d\t'%reaction_export[i][3])
            export_readme.write('Update Number: %-2d\tDetails: '%reaction_export[i][3])
            for j in range(self.__upd_size__):
                export_file.write('%-2d\t%-2d\t'%(reaction_export[i][4][j][0],reaction_export[i][4][j][1]))
                if self.reaction[i][4][j][0]!=0:
                    export_readme.write('[%-8s, %-2d]\t'%(self.reactant[self.reaction[i][4][j][0]][0],reaction_export[i][4][j][1]))
            export_file.write('\n')
            export_readme.write('\n')
        export_file.close()
        export_readme.close()
    
    def save(self, modelName):
        if os.path.exists("./"+modelName)==0:
            os.mkdir(modelName+'_dev')
        save_file=open(modelName+'_dev/reactant','w')
        pickle.dump(self.reactant,save_file)
        save_file.close()
        save_file=open(modelName+'_dev/reaction','w')
        pickle.dump(self.reaction,save_file)
        save_file.close()
        
    def load(self, modelName):
    # no need to add '_dev', or rather, DO NOT add it.
        load_file=open(modelName+'_dev/reactant','r')
        self.reactant=pickle.load(load_file)
        load_file.close()
        load_file=open(modelName+'_dev/reaction','r')
        self.reaction=pickle.load(load_file)
        load_file.close()
        
