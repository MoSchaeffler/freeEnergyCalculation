#!/usr/bin/python3

import numpy as np
import sys
import warnings
import copy
import networkx as nx

def fxn():
    warnings.warn("runtime", RuntimeWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()


class getEnergy:

    def __init__(self,equiProb,probMatrix,nmin,
            temperature=300,timescale=100e-12,
            mts=False,dic=False):
    

        self.nmin = nmin
        self.equiProb = equiProb
        self.probMatrix = probMatrix

        # timescale of the process under study, here time in between frames of the simulation
        self.timescale = timescale 
        #timescale = 1

        # temperature in K
        self.T = temperature
        
        # how to construct ts.data
        # True: For ts.data take mean of forward/backward transition state energy and drop transition states for which only one of both exists.
        # False: use ts1.data and if a transitionstate exists in ts2.data that does not exist in ts1.data use that one
        self.mts = mts
        self.dic = dic

        # run
        self.run()
        
    #############################################
    ####              RUN                   #####
    #############################################
    def run(self):
    
        # load data
        prob = self.loadEqui(self.equiProb)
        trans = self.loadMatrix(self.probMatrix)

        if len(trans) == len(prob):
            nmin = len(trans)
        else:
            print("ERROR: eq. probabilities {} and transition matrix {} of different dimension".format(len(prob),len(trans)))
            sys.exit()
            
        ### Check Data
        self.checkData(prob,trans,self.nmin)
        
        ### calc rate Matrix
        K = self.rateMatrix(trans,self.nmin,self.timescale)
        
        ### calc free energy and transition states
        rhomin, rhots1, rhots2 = self.freeEenergy(K,prob,self.T,self.nmin)
        
        ### write data
        self.write(rhomin, rhots1, rhots2,self.nmin, self.mts,self.dic)
        print("---------------------------------------------------------")
            
            
    #############################################
    ####         Checking Data              #####
    #############################################
    def checkData(self,prob,trans,nmin):
    
    
        print("---------------------------------------------------------")
        print("Checking that equilibrium probabilies are properly normalized")
        s = np.sum(prob)
        
        print("Total Prbability: {}".format(s))
        print("Deviation from 1: {}".format(abs(s-1)))
        
        print("---------------------------------------------------------")
        print("Checking that the maximum probability isn't larger than 1")
        """
        for i in range(nmin):
            for j in range(nmin):
                if trans[i,j] > 1:
                    print(trans[i,j],i,j)
        """
        maxP = np.max(trans)
        print("Maximum probability: {}".format(maxP))

        print("Checking that the sum of all probabilites of movement out of states, are 1")
        sumCheck = np.zeros(nmin)
        for i in range(nmin):
            sumCheck[i] = np.sum(trans[i])

        # calc maximum diviation from 1
        maxDiv = np.max(sumCheck-1)
        print("Max deviation from 1: {:.15f}".format(maxDiv))
        print("---------------------------------------------------------")
        

    #############################################
    ####         Rate Matrix                #####
    #############################################
    def rateMatrix(self,trans,nmin,timescale):

        # M=K-D
        # also change convention for enty k[i,j] from i->j to i<-j to comply with Wales et al.
        # kmatrix = (np.transpose(trans)-np.identity(nmin)) / timescale #(apply timescale later)

        kmatrix = (np.transpose(trans)-np.identity(nmin))/ timescale

        print("Check columns sum to zero")
        check = np.zeros(nmin)
        for i in range(nmin):
            s = np.sum(kmatrix[:,i])
            check[i]=s
            
        print("Max deviation from 0: {:.15f}".format(np.max(np.abs(check))))
        print("---------------------------------------------------------")
        
        return kmatrix    
    
    #############################################
    ####          free energy               #####
    #############################################    
    def freeEenergy(self,kmatrix,prob,T,nmin):
        ### MINIMA in units of kT
        rhomin = -np.log(prob)
        #print(rho)

        kToverh = 1.380649e-23 * T / 6.62607015e-34
        #print("{:e}".format(kToverh))

        ### TRANSITION STATES in units of kT
        rhots1 = np.zeros((nmin,nmin))
        rhots2 = np.zeros((nmin,nmin))
        
        # ignore Runtime Warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fxn()

            for j2 in range(nmin):
                for j1 in range(nmin):
                    rhots2[j1,j2] =  rhomin[j2] - np.log(kmatrix[j1,j2]) + np.log(kToverh)

            for j2 in range(nmin):
                for j1 in range(nmin):
                    rhots1[j1,j2] =  rhomin[j1] - np.log(kmatrix[j2,j1]) + np.log(kToverh)
                    
        return rhomin, rhots1, rhots2


    def effectiveBarrier(self,kmatrix,prob,dic,T,timescale):

        ### symmetrize using average edge capacity
        nmin = len(kmatrix)
        capacity = np.zeros(np.shape(kmatrix)) 

        for i in range(nmin):
            for j in range(i+1,nmin):

                if kmatrix[i,j] != 0 and kmatrix[j,i] != 0:

                    c = (kmatrix[i,j]+kmatrix[j,i])/2
                    capacity[i,j] = c
                    capacity[j,i] = c
                    print(c)

            capacity[i,i] = kmatrix[i,i]
                    
        ### find disconnected states and discard them
        """
        construct the matrix with ts and test if population travels
        trhough network starting in one state by multiplying ts matrix 
        multiple times with initial population
        this samples qualitatively how an initial population probability spreads
        through the network 
        """
                
        # copy capacity into ts matrix
        ts = copy.deepcopy(capacity)
        
        # condidere spread within 15 timesteps
        TS = copy.deepcopy(ts)
        mult = 15
        for i in range(mult):
            TS = np.matmul(TS,ts)
        
        # check spread for initialized population for at each state    
        count_dc = np.zeros(nmin)
        for n in range(nmin):
            v = np.zeros(nmin)
            v[n] = 1
            
            w = np.matmul(TS,v)
            
            nz1 = np.count_nonzero(ts[n])
            nz2 = np.count_nonzero(w)
            
            # a state is onsidered disconnected if fewer states than 
            # the number of multiplications were visisited
            if nz2 <= mult:
                count_dc[n] = 1
                
            #print("State {} has {} connections to other states".format(n,nz1))
            #print("State {} spread to {} states in total".format(n,nz2))
        
        print("Disconnected states: {}".format(np.sum(count_dc)))
        
        if dic == False:
            print("No state dictionary provided")
            print("Mapping of new states to state representation not possible")
            mapstate = False
        else:
            try:
                stateList = self.LoadDict(dic)
                mapstate = True
            except:
                print("Can't read provided dictionary file")
                print(dic)
                mapstate = False
        
        print("Remove disconnected state from network")
        capacity_dc = copy.deepcopy(capacity)
        prob_dc = copy.deepcopy(prob)
        
        
        i = 0 # count loop position of shortened array in respect to loop
        delEntries = []
        stateList_new = []
        stateList_removed = []
        for n in range(nmin):
        
            if mapstate == True:
                state = stateList[n]
            
            # remove dc states
            if count_dc[n] == 1:
                capacity_dc = np.delete(capacity_dc,i,0)
                capacity_dc = np.delete(capacity_dc,i,1)
                prob_dc = np.delete(prob_dc,i)
                delEntries.append(n)
                if mapstate == True:
                    stateList_removed.append([n,state[0]])
                
            else:
                if mapstate == True:
                    stateList_new.append([state[0],i])
            
                i += 1

        ### write new dictionary
        if mapstate == True:    
            with open("mapping.data","w") as f:
                for i in range(len(stateList_new)):
                    f.write("{}\t {}\t \n".format(stateList_new[i][0],stateList_new[i][1]))
                    
            with open("removed.data","w") as f:
                for i in range(len(stateList_removed)):
                    f.write("{}\t {}\t \n".format(stateList_removed[i][0],stateList_removed[i][1]))


        ### build network from new capacity matrix
        G = nx.from_numpy_matrix(capacity_dc)
        print(list(G.nodes))
        ### add the capacity attribute
        for i in range(len(capacity_dc)):
            for j in range(i+1,len(capacity_dc)):
                print(i,j,capacity_dc[i,j])
                G.edges[i, j]['capacity'] = capacity_dc[i,j]
                G.edges[j, i]['capacity'] = capacity_dc[i,j]
        
        # calc min cut using gomory hu
        
        T = nx.gomory_hu_tree(G)
                
        def minimum_edge_weight_in_shortest_path(T, u, v):

            path = nx.shortest_path(T, u, v, weight="weight")

            return min((T[u][v]["weight"], (u, v)) for (u, v) in zip(path, path[1:]))
        
        cut_matrix = np.zeros(np.shape(kmatrix)) 

        for i in range(nmin):
            for j in range(i+1,nmin):

                cut_value, edge = minimum_edge_weight_in_shortest_path(T, i, j)
                print(i,j,cut_value)
                cut_matrix[i,j] = cut_value
                
        return 0
    
    
    

    #############################################
    ####         Helper Functions           #####
    #############################################

    def loadEqui(self,inname):
        with open(inname, 'r') as f:
            lines = f.readlines()
            
            prob = np.array([float(i) for i in lines[0].split()])
        return prob
        
    def loadMatrix(self,inname):
        with open(inname, 'r') as f:
            lines = f.readlines()
            M = []
            for line in lines:
                l = line.split()
                nl = np.array([float(i) for i in l])
                M.append(nl)
            
            
        return np.array(M)
        
    def LoadDict(self,inname):
        dictionary = []
        with open(inname, 'r') as f:
            for line in f:
                if not line.startswith('@'):
                    l = line.split('\t')
                    dictionary.append([eval(l[0]), int(l[1])])

        return dictionary        

    #############################################
    ####          Write Output              #####
    #############################################
    def write(self,rhomin, rhots1, rhots2,nmin,mts,dic=None):
        # minima
        # here the second column is the vibrational contribution 2 ln(2 pi k T/h)
        kb = 1.380649e-23 # J/K
        h = 6.62607015e-34 * 1e12 # J*ps
        vibration = 2*np.log(2*np.pi*kb*self.T/h)
        
        with open("min.data","w") as f:
            for i in range(nmin):
                f.write("{} {}   1  1.0  1.0  1.0 \n".format(rhomin[i],vibration))

        # transitions states
        with open("ts1.data","w") as f:
            for j2 in range(nmin):
                for j1 in range(j2+1,nmin):
                    f.write("{}  0.0  1 {}  {} 1.0  1.0  1.0 \n".format(rhots1[j1,j2],j1+1,j2+1))
                    
        with open("ts2.data","w") as f:
            for j2 in range(nmin):
                for j1 in range(j2+1,nmin):
                    f.write("{}  0.0  1 {}  {} 1.0  1.0  1.0 \n".format(rhots2[j1,j2],j1+1,j2+1))


        if mts == False:
            print("Writing ts.data, using ts1.data and filling up ts with ts2.data")                    
            # merged ts1 amd ts2, dropping inf high transitions
            with open("ts1.data",'r') as f:
                lines1 = f.readlines()

            with open("ts2.data",'r') as f:
                lines2 = f.readlines()

            N = len(lines1)

            rep2 = 0
            with open("ts.data",'w') as f:

                for i in range(N):

                    l1 = lines1[i].replace('"',"")
                    e1 = l1.split()

                    l2 = lines2[i].replace('"',"")
                    e2 = l2.split()

                    if e1[0] == 'inf':

                        if e2[0] != 'inf':

                            f.write(l2)
                            #print(i,l2)
                            rep2 += 1


                    else:
                        f.write(l1)
                        
            print("{} transition states from ts2.data were used".format(rep2))

        if mts == True:
            print("Writing ts.data using mean transition state sheme")
            
            # count discarded transition states and write new ts matrix
            dis1 = 0
            dis2 = 0
            rhots = np.zeros(rhots1.shape)
            for j2 in range(nmin):
                for j1 in range(j2+1,nmin):
                    rho1 = rhots1[j1,j2]
                    rho2 = rhots2[j1,j2]
                    
                    rho = (rho1 + rho2) / 2
                    
                    # count 
                    if rho1 == np.inf and rho2 != np.inf:
                        dis2 += 1
                    if rho1 != np.inf and rho2 == np.inf:
                        dis1 += 1
                        
                    rhots[j1,j2] = rho
                    rhots[j2,j1] = rho
                        
            print("Discarded transition states of ts1.data: {}".format(dis1))
            print("Discarded transition states of ts2.data: {}".format(dis2))
            print("---------------------------------------------------------")
            
            print("Find disconnected states and eliminate from Network")
            
            ### find disconnected states and discard them
            """
            construct the matrix with ts and test if population travels
            trhough network starting in one state by multiplying ts matrix 
            multiple times with initial population
            this samples qualitatively how an initial population probability spreads
            through the network 
            """
                    
            # set ts with non existent ts to 0
            ts = copy.deepcopy(rhots)
            mask = ts == np.inf
            ts[mask] = 0
            
            # condidere spread within 15 timesteps
            TS = copy.deepcopy(ts)
            mult = 15
            for i in range(mult):
                TS = np.matmul(TS,ts)
            
            # check spread for initialized population for at each state    
            count_dc = np.zeros(nmin)
            for n in range(nmin):
                v = np.zeros(nmin)
                v[n] = 1
                
                w = np.matmul(TS,v)
                
                nz1 = np.count_nonzero(ts[n])
                nz2 = np.count_nonzero(w)
                
                if nz2 <= mult:
                    count_dc[n] = 1
                    
                #print("State {} has {} connections to other states".format(n,nz1))
                #print("State {} spread to {} states in total".format(n,nz2))
            
            print("Disconnected states: {}".format(np.sum(count_dc)))
            
            if dic == False:
                print("No state dictionary provided")
                print("Mapping of new states to state representation not possible")
                mapstate = False
            else:
                try:
                    stateList = self.LoadDict(dic)
                    mapstate = True
                except:
                    print("Can't read provided dictionary file")
                    print(dic)
                    mapstate = False
            
            print("Remove disconnected state from network")
            ts_dc = copy.deepcopy(rhots)
            rhomin_dc = copy.deepcopy(rhomin)
            
            
            i = 0 # count loop position of shortened array in respect to loop
            delEntries = []
            stateList_new = []
            stateList_removed = []
            for n in range(nmin):
            
                if mapstate == True:
                    state = stateList[n]
                
                # remove dc states
                if count_dc[n] == 1:
                    ts_dc = np.delete(ts_dc,i,0)
                    ts_dc = np.delete(ts_dc,i,1)
                    rhomin_dc = np.delete(rhomin_dc,i)
                    delEntries.append(n)
                    if mapstate == True:
                        stateList_removed.append([n,state[0]])
                    
                else:
                    if mapstate == True:
                        stateList_new.append([state[0],i])
                
                    i += 1
                
            ### Write new network
            
            nmin = len(rhomin_dc)            
            with open("ts.data","w") as f:
                for j2 in range(nmin):
                    for j1 in range(j2+1,nmin):
                        f.write("{}  0.0  1 {}  {} 1.0  1.0  1.0 \n".format(ts_dc[j1,j2],j1+1,j2+1))
                        
            with open("min.data","w") as f:
                for i in range(nmin):
                    f.write("{} {}   1  1.0  1.0  1.0 \n".format(rhomin_dc[i],vibration))
                
            if mapstate == True:    
                with open("mapping.data","w") as f:
                    for i in range(len(stateList_new)):
                        f.write("{}\t {}\t \n".format(stateList_new[i][0],stateList_new[i][1]))
                        
                with open("removed.data","w") as f:
                    for i in range(len(stateList_removed)):
                        f.write("{}\t {}\t \n".format(stateList_removed[i][0],stateList_removed[i][1]))
            
            






























