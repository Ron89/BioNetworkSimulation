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
	"""
	This module is built to systematically create/modify biochemical networks.
	"""
	__rate_size__=4
	__depend_size__=3
	__upd_size__=3
	__index_length__=6
	reactant={}
	reaction={}
	reactant_const=[]
	reaction_const=[]
	
	def __init__(self):
		random.seed(time.gmtime())
	
	def insert_new_reactant(self,reactant_alias,initial_value=None, unexpected=0):
		"""
		add a new reactant to the network.
		"""
		if unexpected ==1:
			print "reactant",reactant_alias,"not exists, creating instance."
		if initial_value is None: 
			self.reactant[hash(reactant_alias)]=[reactant_alias,int(raw_input("type in initial condition:"))]
		else:
			self.reactant[hash(reactant_alias)]=initial_value
		# (above) change int into float for concentration description
		if unexpected ==1:
			return hash(reactant_alias)
	
	def get_reactant_hash(self,reactant_alias):
		"""
		If reactant already exists in the network, return its hash value as 
		index; if not, trigger insert_new_reactant(...) method to insert it.
		"""
		reactant_hash=hash(reactant_alias)
		if reactant_hash in self.reactant:
			return reactant_hash
		else:
			return self.insert_new_reactant(reactant_alias,unexpected=1)
	
	def list_reactant(self):
		"""list all reactants in the network"""
		print "format:\n index, [reactant name, initial condition]"
		for i in self.reactant:
			print i, self.reactant[i]
	
	def insert_new_reaction(self,code_alias, rate_alias, depend_alias, upd_no_alias, upd_alias):
		"""
		insert a new reaction to the network. Format is as described in README file.
		"""
		# generate index for current reaction
		reaction_index=random.randint(10**self.__index_length__,10**(self.__index_length__+1)-1)
		while reaction_index in self.reaction:
			reaction_index=random.randint(10**self.__index_length__,10**(self.__index_length__+1)-1)
		# generate temporary structure to store reaction
		temp_react=[code_alias,np.zeros(self.__rate_size__),np.zeros(self.__depend_size__, dtype=int),upd_no_alias,[]]
		for i in range(self.__upd_size__):
			temp_react[-1].append([0,0])
		if upd_no_alias<=self.__upd_size__:
			for i in range(upd_no_alias):
				temp_react[4][i][0]=self.get_reactant_hash(upd_alias[i][0])
				temp_react[4][i][1]=upd_alias[i][1]
		else:
			print "exceed maximum updated reactants, fatal"
			return
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
		elif code_alias==3:	 # order 1 promoting Hill equa
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
		elif code_alias==5:   # order 1 inhibitive Hill equa.
			temp_react[1][0] = rate_alias[0]	# k
			temp_react[1][1] = rate_alias[1]	# K
			temp_react[1][2] = rate_alias[2]	# h
			temp_react[2][0] = self.get_reactant_hash(depend_alias[0]) # P
			temp_react[2][1] = self.get_reactant_hash(depend_alias[1]) # factor
		elif code_alias==6:   # inhibitive Hill equa. on MM kinetics
			temp_react[1][0] = rate_alias[0]	# km
			temp_react[1][1] = rate_alias[1]	# Km
			temp_react[1][2] = rate_alias[2]	# Kh
			temp_react[1][3] = rate_alias[3]	# h
			temp_react[2][0] = self.get_reactant_hash(depend_alias[0]) # P
			temp_react[2][1] = self.get_reactant_hash(depend_alias[1]) # factor
			temp_react[2][2] = self.get_reactant_hash(depend_alias[2]) # substrate
		elif code_alias==7:   # simplified MM kinetics
			temp_react[1][0] = rate_alias[0]	# km
			temp_react[1][1] = rate_alias[1]	# Km
			temp_react[2][0] = self.get_reactant_hash(depend_alias[0]) # substrate
		# extend cases here
		else:
			print "what?!"
			print "unknown reaction type, fatal"
			return
		self.reaction[reaction_index]=copy.deepcopy(temp_react)
		
	def list_reaction(self,list_target=None):
		"""
		list all/part of the reactions in the network
		"""
		if list_target==None:
			list_target=self.reaction

		print "format:\n index: reaction type\t reaction rate \n reaction\t reaction rate constants\n"
		for i in list_target:
			# generate RHS of equation
			RHS_reactant={}
			for j in self.reactant:
				RHS_reactant[j]=0
			print i,':',
			if self.reaction[i][0]==0:
				print "order-0\t",
				print "k =", self.reaction[i][1][0]
				print "\t Nul. ->",
			elif self.reaction[i][0]==1:
				print "order-1\t",
				print "k =", self.reaction[i][1][0]
				print "\t",self.reactant[self.reaction[i][2][0]][0],"->",
				RHS_reactant[self.reaction[i][2][0]]+=1
			elif self.reaction[i][0]==2:
				print "order-2\t",
				print "k =", self.reaction[i][1][0]
				print "\t",self.reactant[self.reaction[i][2][0]][0],'+',self.reactant[self.reaction[i][2][1]][0],"->",
				RHS_reactant[self.reaction[i][2][0]]+=1
				RHS_reactant[self.reaction[i][2][1]]+=1
			elif self.reaction[i][0]==3:
				print "Promotive Hill equation\t",
				print "k =", self.reaction[i][1][0], "\tK =", self.reaction[i][1][1], "\tH =", self.reaction[i][1][2]
				print "\t",self.reactant[self.reaction[i][2][0]][0],'+',self.reaction[i][1][2],self.reactant[self.reaction[i][2][1]][0],"<=> [complex] ->",
				RHS_reactant[self.reaction[i][2][0]]+=1
				RHS_reactant[self.reaction[i][2][1]]+=self.reaction[i][1][2]
			elif self.reaction[i][0]==4:
				print "MM kinetics\t",
				print "k =", self.reaction[i][1][0], "\tK =", self.reaction[i][1][1]
				print "\t",self.reactant[self.reaction[i][2][0]][0],'+',self.reactant[self.reaction[i][2][1]][0],"<=> [complex] ->",
				RHS_reactant[self.reaction[i][2][0]]+=1
				RHS_reactant[self.reaction[i][2][1]]+=1
			elif self.reaction[i][0]==5:
				print "Inhibitive Hill equation\t",
				print "k =", self.reaction[i][1][0], "\tK =", self.reaction[i][1][1], "\tH =", self.reaction[i][1][2]
				print "\t","[complex] <=>",self.reactant[self.reaction[i][2][0]][0],'+',self.reaction[i][1][2],self.reactant[self.reaction[i][2][1]][0], "->"
				RHS_reactant[self.reaction[i][2][0]]+=1
				RHS_reactant[self.reaction[i][2][1]]+=self.reaction[i][1][2]
			elif self.reaction[i][0]==6:
				print "Inhibitive Hill equation on MM kinetics\t",
				print "km =", self.reaction[i][1][0], "\tKm =", self.reaction[i][1][1], "\tKh =", self.reaction[i][1][2], "\tH =", self.reaction[i][1][3]
				print "\t", self.reactant[self.reaction[i][2][2]][0],"+ [complex] <=>",self.reactant[self.reaction[i][2][2]][0],"+",self.reactant[self.reaction[i][2][0]][0],'+',self.reaction[i][1][3],self.reactant[self.reaction[i][2][1]][0], "->",
				RHS_reactant[self.reaction[i][2][0]]+=1
				RHS_reactant[self.reaction[i][2][1]]+=self.reaction[i][1][3]
				RHS_reactant[self.reaction[i][2][2]]+=1
			elif self.reaction[i][0]==7:
				print "Simplified MM kinetics\t",
				print "v = k[E] =", self.reaction[i][1][0], "\tK =", self.reaction[i][1][1]
				print "\t",self.reactant[self.reaction[i][2][0]][0],'+ [E] <=> [complex] -> [E] +',
				RHS_reactant[self.reaction[i][2][0]]+=1
			# extend cases here
			else:
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
		"""
		export the network into a format(folder with reactant, reaction, readme 
		file in it) recognizable by coarse grained simulation modules defined in 
		coarseGrainedOperation.h.
		"""
		if os.path.exists("./"+exportName)==0:
			os.mkdir(exportName)
		reactant_export=copy.deepcopy(self.reactant)
		reaction_export=copy.deepcopy(self.reaction)
		index_temp=0
		export_readme=open(exportName+'/readme','w')
		export_readme.write("Exported time: %2d/%2d/%4d.\n\n"%(time.localtime().tm_mday,time.localtime().tm_mon,time.localtime().tm_year))
		export_readme.write("Reactants:\n")
		export_file=open(exportName+'/reactant','w')
		for i in sorted(reactant_export):
			reactant_export[i].append(index_temp)
			index_temp+=1
			export_file.write('%d\n'%reactant_export[i][1])
			export_readme.write('%-2d\t%-8s\t%d\n'%(reactant_export[i][-1],reactant_export[i][0],reactant_export[i][1]))
		export_file.close()
		export_readme.write("\nReactions\n")
		export_file=open(exportName+'/reaction','w')
		index_temp=0
		for i in sorted(reaction_export):
			if reaction_export[i][0]<=5:
				for j in range(min(reaction_export[i][0],2)):
					reaction_export[i][2][j]=reactant_export[reaction_export[i][2][j]][-1]
			elif reaction_export[i][0]==6:
				for j in range(3):
					reaction_export[i][2][j]=reactant_export[reaction_export[i][2][j]][-1]
			elif reaction_export[i][0]==7:
				for j in range(1):
					reaction_export[i][2][j]=reactant_export[reaction_export[i][2][j]][-1]
			for j in range(reaction_export[i][3]):
				reaction_export[i][4][j][0]=reactant_export[reaction_export[i][4][j][0]][-1]
			export_readme.write("%-4d\t"%(index_temp))
			index_temp+=1
			export_file.write('%1d\t'%reaction_export[i][0])
			export_readme.write('Type: %1d\t'%reaction_export[i][0])
			for j in range(self.__rate_size__):
				export_file.write('%.4e\t'%reaction_export[i][1][j])
#				export_readme.write('%.4e\t'%reaction_export[i][1][j])
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
		"""
		save the current network into a folder. Folder name will be 
		modelName+'_dev'. pickle is required for the function.
		"""
		if os.path.exists("./"+modelName+'_dev')==0:
			os.mkdir(modelName+'_dev')
		save_file=open(modelName+'_dev/reactant','w')
		pickle.dump(self.reactant,save_file)
		save_file.close()
		save_file=open(modelName+'_dev/reaction','w')
		pickle.dump(self.reaction,save_file)
		save_file.close()
		save_file=open(modelName+'_dev/reaction_const','w')
		pickle.dump(self.reaction_const,save_file)
		save_file.close()
		save_file=open(modelName+'_dev/reactant_const','w')
		pickle.dump(self.reactant_const,save_file)
		save_file.close()
		
	def load(self, modelName):
		""" 
		load a model from previously saved files. Only model name is required to 
		load it. No need to add '_dev', or rather, DO NOT add it.
		"""
		load_file=open(modelName+'_dev/reactant','r')
		self.reactant=pickle.load(load_file)
		load_file.close()
		load_file=open(modelName+'_dev/reaction','r')
		self.reaction=pickle.load(load_file)
		load_file.close()
	# constant marker not exist in initial version, check if available before start reading
		if os.path.isfile(modelName+'_dev/reaction_const'):
			load_file=open(modelName+'_dev/reaction_const','r')
			self.reaction_const=pickle.load(load_file)
			load_file.close()
		else:
			self.const_reaction={}
		if os.path.isfile(modelName+'_dev/reactant_const'):
			load_file=open(modelName+'_dev/reactant_const','r')
			self.reactant_const=pickle.load(load_file)
			load_file.close()
		else:
			self.const_reactant={}
	# in case the saved file is from previous versions, translate the structure into current standard
	# * update target: list -> array
	#	for i in self.reaction:
	#		if self.reaction[i][-1]==list:
	#			self.reaction[i][-1]=np.array(self.reaction[i][-1])
	
	def scaling(self, eta_c_alias, eta_t_alias):
		"""
		scale the network by a certain spacial(eta_c_alias) and time(eta_t_alias)
		ratio. reactions and reactants noted as constant will not be subjected to 
		scaling. Non-constant reactions containing constant reactants will adjust 
		their reaction rates accordingly.
		"""
		for i in self.reactant:
			if i not in self.reactant_const:
				self.reactant[i][1]*=eta_c_alias
		for i in self.reaction:
			if i not in self.reaction_const:
			# dimension determination
				rate_dimension=np.zeros((2,self.__rate_size__))  # first row: spatial dimention, second row: time dimension
				rate_dimension[0][0]= (1 if self.reaction[i][0]==0 else (-1 if self.reaction[i][0]==2 else 0))
				rate_dimension[0][1]= 1 if (self.reaction[i][0]==3 or self.reaction[i][0]==4) else 0
				rate_dimension[1][0]= -1;
				# adjust dimension according to constant reactants
				if self.reaction[i][0]==1:
					if self.reaction[i][2][0] in self.reactant_const:
						rate_dimension[0][0]+=1
				if self.reaction[i][0]==2:
					if self.reaction[i][2][0] in self.reactant_const:
						rate_dimension[0][0]+=1
					if self.reaction[i][2][1] in self.reactant_const:
						rate_dimension[0][0]+=1
				if self.reaction[i][0]==3:
					if self.reaction[i][2][0] in self.reactant_const:
						rate_dimension[0][0]+=1
				if self.reaction[i][0]==4:
					if self.reaction[i][2][1] in self.reactant_const:
						rate_dimension[0][0]+=1
				# scale the reaction
				for j in range(self.__rate_size__):
					self.reaction[i][1][j]*=eta_c_alias**rate_dimension[0][j]*eta_t_alias**rate_dimension[1][j]

	def knockout_reactant(self, knockout_list):
		"""
		remove the listed reactants from the network. Reactants in the list need to be hashed value.
		* all manual hashing should be changed in future versions.
		"""
		knockout_reaction_list=set()
		action_confirmed=1;
		for target in knockout_list:
			for affectTest in self.reaction:
				if (target in self.reaction[affectTest][2]) or (target in np.array(self.reaction[affectTest][-1])):
					knockout_reaction_list.add(affectTest)
		if len(knockout_reaction_list)!=0:
			action_confirmed=0
			print "Following reactions are affected:"
			print '\t',knockout_reaction_list
			print "Detailed List:"
			self.list_reaction(knockout_reaction_list)
			print "Reactions listed above will also be removed from the network by the action"
			if raw_input("Type \'yes\' to confirm deletion: ")=="yes":
				action_confirmed=1

		if action_confirmed==1:
			for target in knockout_reaction_list:
				self.reaction.pop(target)
				if target in self.reaction_const:
					self.reaction_const.pop(target)
			for target in knockout_list:
				self.reactant.pop(target)
				if target in self.reactant_const:
					self.reactant_const.pop(target)

	def duplicate_reaction(self, duplicate_target):
		"""
		create a duplicate reaction in the network. return the index(key) of the new duplicated reaction.
		"""
		# generate index for current reaction
		reaction_index=random.randint(10**self.__index_length__,10**(self.__index_length__+1)-1)
		while reaction_index in self.reaction:
			reaction_index=random.randint(10**self.__index_length__,10**(self.__index_length__+1)-1)
		self.reaction[reaction_index]=copy.deepcopy(self.reaction[duplicate_target])
		return reaction_index;
