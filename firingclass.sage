
class FiringClass():
    '''
    A class that generates, holds, and displays all information about a certain L,M firing setup.
    '''
    
    __sorting = 'groups'          #Sorting metric
    __is_cycle = False            #Cycle flag
    __positive = False            #Positive flag
    __computed = False            #Computed flag
    
    G = None                      #Graph object
    L = None                      #Laplacian Matrix of Gphi
    Linv = None                   #Inverse of L
    M = None                      #Laplacian Matrix of |G|
    Minv = None                   #Inverse of M
    LMinv = None                  #LM-1
    MLinv = None                  #ML-1
    n = 0                         #Dimension of G
    dL = 0                        #Determinant of L
    dM = 0                        #Determinant of M
    bucket_size = 0               #Size of each fractional bucket
    cmax = None                   #Maximal Critical Vector
    dmax = None                   #Degree Vector
    sandpile = None               #Sandpile object
    group_structure = None        #List of interesting invariant factors
    
    supers = []                   #Supers for |G|
    supers_signed = []            #Supers in S+
    supers_preimage = []          #Preimages in R+
    supers_frac = []              #Fractional parts of the preimages
    supers_groups = []            #Grouping Lists (with repeats)
    supers_groups_d = []          #Grouping Map (with counts)

    criticals = []                #Crits for |G|
    criticals_signed = []         #Crits in S+
    criticals_preimage = []       #Preimages in R+
    criticals_frac = []           #Fractional parts of the preimage
    criticals_groups = []         #Grouping Lists (with repeats)
    criticals_groups_d = []       #Grouping Map (with counts)

    
    def __init__(self, L, G):
        '''
        The constructor FiringClass(L,G), generates all shallow data of the signed graph.
        Computes LM-1, cmax, sandpiles, etc
        
        Arguments:
        L - The Laplacian matrix of G
        G - A graph object
        '''
        self.G = copy(G)  

        self.L = copy(L)
        
        self.Linv = self.L.inverse()
        self.M = self.G.laplacian_matrix()[:-1,:-1]
        self.Minv = self.M.inverse()

        self.__positive = self.L == self.M
        self.__is_cycle = isCycle(self.M)
        
        self.LMinv = self.L*self.Minv
        self.MLinv = self.M*self.Linv
        
        self.dL = abs(self.L.determinant())
        self.dM = abs(self.M.determinant())
        
        self.n = Integer(len(self.G))        
        self.cmax = get_cmax(self.M)
        self.dmax = vector([c+1 for c in self.cmax])

        lL = (self.Linv)*IntegralLattice(self.n-1)
        lM = (self.Minv)*IntegralLattice(self.n-1)
        both = lM.intersection(lL)
        quo = both.quotient(both.intersection(IntegralLattice(self.n-1)))
        inv = quo.invariants()
        self.bucket_size = list_prod(inv)
        
        self.group_structure = invariantFactors(self.L)
        
        self.sandpile = Sandpile(self.G,self.n)

        self.__computed = False
        self.__sorted = 'groups'
        
        self.supers = []                   
        self.supers_signed = []            
        self.supers_preimage = []     
        self.supers_frac = []
        self.supers_groups = []            
        self.supers_groups_d = []          

        self.criticals = []                
        self.criticals_signed = []         
        self.criticals_preimage = [] 
        self.criticals_frac
        self.criticals_groups = []         
        self.criticals_groups_d = []       

        
    def __cone_integrals_full(self):
        '''
        Computes the integral points in the cone S+ bounded by the degree vector dmax. Unoptimized.
        This is only here so 2025 has an easier time learning how we did it
        '''
        inequalities_bounding_cube = [[self.dmax[i]]+list(-matrix.identity(self.n-1)[i]) for i in range(self.n-1)] + [[0]+list(matrix.identity(self.n-1)[i]) for i in range(self.n-1)]
    
        prepoly = Polyhedron(ieqs=inequalities_bounding_cube)
        newpoly = self.LMinv*prepoly
        cone = Cone(self.LMinv.transpose()).intersection(newpoly)

        cone_integrals = cone.integral_points()
        return cone_integrals
        
        
    def __cone_integrals_split(self):
        '''
        Computes the integral points in the cone S+ bounded by the degree vector dmax.
        Uses a splitting optimization that reduces the dimension of the polyhedra by the number of positive incident vertices.
        Works best when the dimension is reduced by about 50%
        '''
        positiveVerts = []
        eye = identity_matrix(self.n-1)
        
        for i in range(self.n -1):
            if self.LMinv[i] == eye[i]:
                positiveVerts += [i+1]
                
        if not positiveVerts or len(positiveVerts) == self.n-1:
            #optimization does nothing on all negative vertices and all positive graphs
            return self.__cone_integrals_full()

        #Small optimization: use the vertices with the smallest degree so the number of spaces scales slower. 
        positiveVerts.sort(key = lambda x: self.dmax[x-1])      
        positiveVerts = positiveVerts[:(self.n-1)//2] #small optimization: use around half of the positive vertices to keep balance between exponentials. ~2x speedup
        positiveVerts.sort()
        
                
        bound = vector([ self.dmax[i] for i in range(len(self.dmax)) if not i+1 in positiveVerts])
        Vertices = [bound.pairwise_product(vector(mask)) for mask in product([0,1], repeat = len(bound))]
        
        for i in range(len(Vertices)):
            Vertices[i] = list(Vertices[i])
            for j in positiveVerts:
                Vertices[i].insert(j-1,0)
            Vertices[i] = vector(Vertices[i])
            
        baseP = self.LMinv*Polyhedron(Vertices)
        
        offsets = [zero_vector(self.n-1)]
        
        for dim in positiveVerts:
            newOffs = []
            for off in offsets:
                for c in range(1, self.dmax[dim-1]):
                    newOffs.append(off + matrix.identity(self.n-1)[dim-1] * c)
            offsets += newOffs
                 
        pList = [baseP + self.LMinv*offset for offset in offsets]
        
        cone = Cone(self.LMinv.transpose())

        cone_integrals = []
        @parallel(24)
        def get_int_list(p):
            return (cone&p).integral_points()
        
        cone_integrals2 = get_int_list(pList)
        cone_integrals = []
        for a in cone_integrals2:
            cone_integrals += list(a[1]) 
        return cone_integrals

    def compute_supers_p_only(self, supers, method = 'split'):
        '''
        Computes only the preimage superstables of the graph. Used in the vertex switching optimization.
        
        Arguments:
        supers -- The superstable configurations of the unsigned graph
        method -- The method to be used to construct the polyhedron. 'full' or 'split'
        '''
        
        if self.__positive:
            self.supers_preimage = supers
            return
        
        #Finding integral points
        
        upper_bound_rat = self.dmax - (1/self.dL)*vector([1]*(self.n - 1))

        if method == 'full':
            cone_integrals = self.__cone_integrals_full()
        else:
            cone_integrals = self.__cone_integrals_split()
        
        cone_preimages = []

        
        #Filtering integral points
        for v in cone_integrals:
            preimg = self.MLinv*v
            fv = floor_vec(preimg)
            if fv in supers:
                self.supers_preimage.append(preimg)
        
        
    def compute_sss(self,method = 'split', t = False, new = False, sort = True):
        '''
        Computes the supercritical configurations of the graph, along with the group structure.
        Finds the switching equivalent graph that is fastest to compute, and computes it.
        
        Arguments:
        method -- The method to be used to construct the polyhedra.
                  'full' creates one n dimensional polyhedra and checks the entire thing. Used for negative graphs
                  'split' splits the polyhedron based on positive vertices. Works better on graphs with few negative edges. 
                   automatically switches to full if the optimization will not help (i.e. no positive vertices)
        t      -- Displays the time taken to compute each step and the total amount of time taken. 
        '''
        self.__computed = True
        
        if t:
            timeT = time.time()
            timeTotal = time.time()
        
        
        #Unsigned Supercriticals
        if self.__is_cycle:
            self.criticals = [vector(s) for s in Permutations([0] + [1]*(self.n-1), self.n-1)]
            self.supers = [equivalentSuperstable(self.sandpile, s) for s in self.criticals]
        else:
            self.criticals = self.sandpile.recurrents(False)
            self.supers = [equivalentSuperstable(self.sandpile,c) for c in self.criticals]
            self.criticals = [vector(s) for s in self.criticals]
        
        if t:
            timeT = timeSince(timeT, True, "SSS - Original supers")
        
        if self.__positive:
            self.supers_signed = self.supers
            self.supers_preimage = self.supers
            self.supers_groups = self.supers
            self.supers_groups_d = self.supers
            self.supers_frac = [frac_vec(v) for v in self.supers_preimage]
        
            self.criticals_signed = self.criticals
            self.criticals_preimage = self.criticals
            self.criticals_groups = self.criticals
            self.criticals_groupds_d = self.criticals
            self.criticals_frac = [frac_vec(v) for v in self.criticals_preimage]

            self.sort(self.__sorting)
            return
        
        #sign switching optimization
        Li,path = None, None
        
        if self.__is_cycle:
            Li = self.M
            path = getCycleSwitches(self.L,self.M)
        elif self.dL == self.dM:
            Li = self.M
            path = getSwitches(self.L, self.M, G)
        else:
            Li, path = getSwitchRep(self.L,self.G, 'posVerts')  
        
        
        if t:
            timeT = timeSince(timeT, True, "SSS - Get switch rep")
        
        swapped_preimages = []
        
        if Li == self.M:
            swapped_preimages = self.supers
        else:
            fc = FiringClass(Li,G)
            fc.compute_supers_p_only(self.supers, method)
            swapped_preimages = fc.supers_preimage
            
        
        if t:
            timeT = timeSince(timeT, True, "SSS - Compute switched supers")
        
        E = matrix.identity(self.n-1, self.n-1)
        for a in path:
            E[a-1,a-1] = -1
        delta = self.M*E*self.Minv


        swapped_preimages = [ delta*v for v in swapped_preimages]
        swapped_preimages = [ (vector(get_sp_config(self.sandpile, floor_vec(v)).equivalent_recurrent().values()), frac_vec(v)) for v in swapped_preimages]
            
        superindexmap = {tuple(s):i for i, s in enumerate(self.criticals)}
        self.criticals_preimage = [v + s for s, v in swapped_preimages]
        self.supers_preimage = [ v + self.supers[superindexmap[tuple(s)]] for s, v in swapped_preimages]
        
        self.supers_signed = [self.LMinv * v for v in self.supers_preimage]
        self.criticals_signed = [self.LMinv * v for v in self.criticals_preimage]
        
        if t:
            timeT = timeSince(timeT, True, "SSS - Switch supers back")
        
        #Supercritical Groupings
        self.supers_groups = [floor_vec(v) for v in self.supers_preimage]
        self.criticals_groups = [floor_vec(v) for v in self.criticals_preimage]
        
        ss_groups = count_sublists(self.supers_groups)
        crit_groups = count_sublists(self.criticals_groups)
        self.supers_groups_d = [(k, ss_groups[k]) for k in ss_groups]
        self.criticals_groups_d = [(k, crit_groups[k]) for k in crit_groups]
        
        self.supers_frac = [frac_vec(v) for v in self.supers_preimage]
        self.criticals_frac = [frac_vec(v) for v in self.criticals_preimage]
        
        if t:
            timeSince(timeTotal, label = "SSS - Total Runtime")
            print()
            
        if sort:
            self.sort(self.__sorting)
        
    def displayGraph(self):
        #Visually displays the graph
        displayGraph(self.L,self.G)
        
        
    def displayShallow(self):
        #Displays shallow information about the graph such as the group structure
        self.displayGraph()
        print(f"Cardinality of unsigned critical group: {self.dM}")
        print(f"Cardinality of signed critical group: {self.dL}")
        print(f"\nNumber of Switching Equivalence Classes: {numSwitchClasses(self.G)}")
        print(f"Size of Each Switching Equivalence Class: {sizeSwitchClasses(self.G)}")
        print(f"\nGroup structure: {getGroupStructure(self.L)}")
        print(f"\nBucket size: {self.bucket_size}")
        print("\nLM-1:")
        print(self.LMinv)
        
        
        
    def displayAll(self, groups = True, limit = True):
        #Displays all information about the graph including the supercriticals in R+ and S+
        self.displayShallow()
        
        if not self.__computed:
            self.compute_sss()
        print(f"\nSorted by {self.__sorting}:")
        print("\nOriginal Superstables and Criticals:")
        columnDisplay([self.supers, self.criticals], ["Super", "Critical"], limit = limit)
        
        if self.__positive:
            return

        if groups:
            print("\nSigned Superstables and Criticals in S+:")
            columnDisplay([self.supers_signed,self.supers_groups, self.criticals_signed, self.criticals_groups], ["Super","Group", "Critical", "Group"], limit = limit)
        
            print("\nPreimages in R+:")
            columnDisplay([self.supers_preimage,self.supers_groups, self.criticals_preimage, self.criticals_groups], ["Super","Group", "Critical","Group"], limit=limit)
        else:
            print("\nSigned Superstables and Criticals in S+:")
            columnDisplay([self.supers_signed, self.criticals_signed], ["Super","Critical"])
        
            print("\nPreimages in R+:")
            columnDisplay([self.supers_preimage, self.criticals_preimage,], ["Super", "Critical"])
        
        print("\nPreimage Groupings:")
        columnDisplay([[a for a, b in self.supers_groups_d],  [a for a,b in self.criticals_groups_d],[b for a,b in self.supers_groups_d], [b for a,b in self.criticals_groups_d]],["Super", "Critical", "scnt", "ccnt"], limit = True)
        
    def displayGroupStructure(self):
        #Displays the group structure of the signed graph
        print(f"\nGroup structure: {getGroupStructure(self.L)}")
        
    def displayNumSwitchClasses(self):
        #Displays the number and size of the switching equivalence classes of G
        print(f"\nNumber of Switching Equivalence Classes: {numSwitchClasses(self.G)}")
        print(f"\nSize of Each Switching Equivalence Class: {sizeSwitchClasses(self.G)}")
        
        
    def displayUnsigned(self, limit=False):
        #Displays the supercriticals of the unsigned graph G
        if not self.__computed:
            self.compute_sss()
        print("\nOriginal Superstables and Criticals:")
        columnDisplay([self.supers, self.criticals], ["Super", "Critical"], limit = limit)
        
    def displaySigned(self, groups = True, limit = False):
        #Displays the supercriticals of the signed graph in S+
        if not self.__computed:
            self.compute_sss()
            
        if groups:
            print("\nSigned Superstables and Criticals in S+:")
            columnDisplay([self.supers_signed,self.supers_groups, self.criticals_signed, self.criticals_groups], ["Super","Group", "Critical", "Group"], limit = limit)
        else:
            print("\nSigned Superstables and Criticals in S+:")
            columnDisplay([self.supers_signed, self.criticals_signed], ["Super", "Critical"], limit = limit)
            
    def displayPreimages(self, groups = True, limit = False):
        #Displays the supercriticals of the signed graph in R+
        if not self.__computed:
            self.compute_sss()
        if groups:
            print("\nPreimages in R+:")
            columnDisplay([self.supers_preimage,self.supers_groups, self.criticals_preimage, self.criticals_groups], ["Super","Group", "Critical","Group"], limit=limit)
        else:
            print("\nPreimages in R+:")
            columnDisplay([self.supers_preimage, self.criticals_preimage], ["Super","Critical"], limit=limit)
            
        
    def displayGroupings(self, limit=False):
        #Displays the groupings of the supercriticals and the number of entries in each group
        if not self.__computed:
            self.compute_sss()
        print("\nPreimage Groupings:")
        columnDisplay([[a for a, b in self.supers_groups_d],  [a for a,b in self.criticals_groups_d],[b for a,b in self.supers_groups_d], [b for a,b in self.criticals_groups_d]],["Super", "Critical", "scnt", "ccnt"], limit = limit)
        
    def displayFractionalParts(self):
        #Displays the fractional parts of the supercriticals in R+
        if not self.__computed:
            self.compute_sss()
        print("\nFractional Parts:")
        columnDisplay([self.supers_frac, self.criticals_frac], ["Supers", "Criticals"])

    def displayBijection(self):
        '''
        Displays the bijection between the superstables and criticals given by the phimap
        '''
        if not self.__computed:
            self.compute_sss()
            
        involution = {} # involution on unsigned superstables

        # finding fixed points
        fixed_points = []
        for s in self.supers:
            if integer_list(self.LMinv*(self.cmax - 2*s)):
                fixed_points += [s]
                involution[tuple(s)] = s

        # finding the rest of the involution
        for s1 in self.supers:
            if tuple(s1) not in involution:
                for s2 in self.supers:
                    if frac_vec(self.Minv*(self.cmax-s1)) == frac_vec(self.Minv*s2):
                        involution[tuple(s1)] = s2
                        involution[tuple(s2)] = s1
                        break
        
        # computing bijection
        computed = []
        for sp in self.supers_preimage:
            computed += [frac_vec(sp) + self.cmax-involution[tuple(floor_vec(sp))]]
        
        # printing results
        print("\nBijection: ")
        columnDisplay([self.supers_preimage, computed], ["Super", "Critical"])

        inv = []
        used = set()
        for k in involution:
            if k not in used:
                inv += [(k,tuple(involution[k]))]
                used.add(k)
                used.add(tuple(involution[k]))
        
        print("\nInvolution Pairs:")
        columnDisplay([inv], [":"])
        print()

    def sort(self, metric = 'groups', index = 0): 
        '''
        Sorts the supercriticals by a metric.
        
        Arguments:
        metric -- The metric which the supercriticals will be sorted by. Options are shown below
        index  -- Optional index to sort by if the metric is in size order.
        '''
        metric = metric.lower()
        self.__sorting = metric
        ssZip = list(zip(self.supers_signed, self.supers_preimage, self.supers_groups, self.supers_frac))
        crZip = list(zip(self.criticals_signed, self.criticals_preimage, self.criticals_groups, self.criticals_frac))
        
        
        if metric == 'preimg':
            #sort by size order of the preimages
            self.supers.sort()
            self.criticals.sort()
            if self.__positive: 
                return
            ssZip.sort(key = lambda x: (x[1][index],x[1]))
            crZip.sort(key = lambda x: (x[1][index],x[1]))
            self.supers_groups_d.sort(key=lambda x: x[0])
            self.criticals_groups_d.sort(key=lambda x:x[0])
        elif metric == 'image':
            #sort by size order of the images
            self.supers.sort()
            self.criticals.sort()
            if self.__positive: 
                return
            ssZip.sort(key = lambda x: (x[0][index],x[0]))
            crZip.sort(key = lambda x: (x[0][index],x[0]))
            self.supers_groups_d.sort(key=lambda x: x[0])
            self.criticals_groups_d.sort(key=lambda x:x[0])
        elif metric == 'frac':
            #sort by size order of the fractional parts
            self.supers.sort()
            self.criticals.sort()
            if self.__positive: 
                return
            ssZip.sort(key = lambda x: (x[3][index],x[3]))
            crZip.sort(key = lambda x: (x[3][index],x[3]))
            self.supers_groups_d.sort(key=lambda x: x[0])
            self.criticals_groups_d.sort(key=lambda x:x[0])
        elif metric == 'duality':
            #sort ss by group size order, and match groups by duality
            self.supers.sort()
            dMap = {tuple(self.cmax - self.supers[i]):i for i in range(len(self.supers))}
            self.criticals.sort(key = lambda x: dMap[tuple(x)])            
            if self.__positive: 
                return
            ssZip.sort(key = lambda x: dMap[tuple(self.cmax - x[2])])
            crZip.sort(key = lambda x: dMap[tuple(x[2])])
            
            self.supers_groups_d.sort(key=lambda x: dMap[tuple(self.cmax - vector(x[0]))])
            self.criticals_groups_d.sort(key=lambda x: dMap[x[0]])
            
        else:
            #sort by size order of groups, matching supers and criticals
            self.supers.sort()
            gMap = {tuple(get_sp_config(self.sandpile, self.supers[i]).equivalent_recurrent().values()):i for i in range(len(self.supers))}
            self.criticals.sort(key = lambda x: gMap[tuple(x)])
            if self.__positive: 
                return
            ssZip.sort(key = lambda x: (x[2],x[1]))
            crZip.sort(key = lambda x: (gMap[tuple(get_sp_config(self.sandpile, x[2]).values())],x[1]))

            self.supers_groups_d.sort(key=lambda x: (x[1], gMap[tuple(get_sp_config(self.sandpile,x[0]).equivalent_recurrent().values())]))
            self.criticals_groups_d.sort(key=lambda x: (x[1], gMap[tuple(get_sp_config(self.sandpile,x[0]).equivalent_recurrent().values())]))
        
        
        
        self.supers_signed = [ss for (ss, ps, gs, fs), (sc, pc, gc, fc) in zip(ssZip, crZip)]
        self.supers_preimage = [ps for (ss, ps, gs, fs), (sc, pc, gc, fc) in zip(ssZip, crZip)]
        self.supers_groups = [gs for (ss, ps, gs, fs), (sc, pc, gc, fc) in zip(ssZip, crZip)]
        self.supers_frac = [fs for (ss, ps, gs, fs), (sc, pc, gc, fc) in zip(ssZip, crZip)]
        self.criticals_signed = [sc for (ss, ps, gs, fs), (sc, pc, gc, fc) in zip(ssZip, crZip)]
        self.criticals_preimage = [pc for (ss, ps, gs, fs), (sc, pc, gc, fc) in zip(ssZip, crZip)]
        self.criticals_groups = [gc for (ss, ps, gs, fs), (sc, pc, gc, fc) in zip(ssZip, crZip)]
        self.criticals_frac = [fc for (ss, ps, gs, fs), (sc, pc, gc, fc) in zip(ssZip, crZip)]
        

    def getBuckets(self):
        '''
        Returns a map from the fractional parts to the supers and criticals that correspond to it.
        
        d[tuple(frac)] = [ [supers], [criticals] ]
        '''
        if not self.__computed:
            self.compute_sss()
            
        sd = {}
        for s in self.supers_preimage:
            sd[tuple(frac_vec(s))] = sd[tuple(frac_vec(s))] + [s] if tuple(frac_vec(s)) in sd else [s]

        cd = {}
        for c in self.criticals_preimage:
            cd[tuple(frac_vec(c))] = cd[tuple(frac_vec(c))] + [c] if tuple(frac_vec(c)) in cd else [c]

        d = {}
        for k in sd:
            d[k] = [sd[k], cd[k]]
        return d
        
    def vertexSwitch(self, vert):
        #Returns the firingclass that corresponds to the vertex switch of this one
        Lnew = copy(self.L)
        Gnew = copy(self.G)
        Lnew = vertexSwitch(Lnew, vert)
        return FiringClass(Lnew,Gnew)
    
    def vertexSwitches(self, V):
        #Returns the firingclass that corresponds to the vertex switches of this one
        Lnew = copy(self.L)
        Gnew = copy(self.G)
        Lnew = vertexSwitches(Lnew, V)
        return FiringClass(Lnew,Gnew)

    def displayPoset(self,figsize = 25):
        #Displays the poset given by elementwise >= on the given list "vectors" in Rn
        # Colors will be s-frackets and compliments, and the labels are printed underneath
        if not self.__computed:
            self.compute_sss(sort = False)
            
        cmax = self.LMinv*self.cmax
        
        edges = []
        images = [self.LMinv*s for s in self.supers]
        fracs = [tuple(frac_vec(s)) for s in images]
        fracsNum = {v: i for i,v in enumerate(list(set(fracs)))} 
        
        fixedColor = {"black":[]}

        from sage.all import colors
        
        colors = list_subtract(list(colors)[1:],['black'])
        colen = len(colors)
        fracColors = {}
        
        for i in range(len(fracs)):
            if fracs[i] not in fracColors:
                fracColors[fracs[i]] = colors[(i+1)%colen]
                fracColors[tuple(frac_vec(cmax-vector(fracs[i])))] = colors[(i+1)%colen]
                fixedColor[colors[(i+1)%colen]] = []

        for s1 in images:
            for s2 in images:
                if le_vec(s1, s2):
                    edges += [[str(self.MLinv*s1) + "\n"+ str(fracsNum[tuple(frac_vec(s1))]) + "\n\n\n\n", str(self.MLinv*s2) + "\n" + str(fracsNum[tuple(frac_vec(s2))]) + "\n\n\n\n"]]
            if frac_vec(2*s1) == frac_vec(cmax):
                fixedColor["black"] += [str(self.MLinv*s1) + "\n"+ str(fracsNum[tuple(frac_vec(s1))]) + "\n\n\n\n"]
            else:
                fixedColor[fracColors[tuple(frac_vec(s1))]] += [str(self.MLinv*s1) + "\n"+ str(fracsNum[tuple(frac_vec(s1))]) + "\n\n\n\n"]

        digraph = DiGraph(edges, loops=True).transitive_reduction()
        poset = Poset(digraph)
        plot = poset.plot(figsize = figsize, edge_color='lightgray', vertex_colors = fixedColor)
        
        done = set()
        frac_pairs = []
        for s1 in images:
            if tuple(frac_vec(s1)) not in done:
                frac_pairs += [(fracsNum[tuple(frac_vec(s1))], fracsNum[tuple(frac_vec(cmax-s1))])]
                done.add(tuple(frac_vec(s1)))
                done.add(tuple(frac_vec(cmax-s1)))
        print("Fracket Indices: ", frac_pairs)
        
        plot.show()
        print(latex(plot))

    
    def isPositive(self):
        #Checks if L == M
        return self.__positive

    def switchingEquivalent(self, other):
        #Checks if two firingclasses are switching equivalent    
        G = copy(self.G)
        L = copy(self.L)
        Lo = copy(other.L)
        n = G.num_verts()
        indexMap = []
        switchMap = [Bitset() for i in range(n)]
    
        for i in range(n-1):
            for j in range(i+1, n-1):
                if L[i,j] != 0:
                    indexMap += [(i,j)]
                    switchMap[i+1].add(len(indexMap)-1)
                    switchMap[j+1].add(len(indexMap)-1)     
        l = len(indexMap)
    
        current = []
        for i in range(l):
            if(L[indexMap[i]] == 1):
                current.append(i)
        current = FrozenBitset(current,capacity = l)
    
        bito = []
        for i in range(l):
            if(Lo[indexMap[i]] == 1):
                bito.append(i)
        bito = FrozenBitset(bito,capacity = l)
    
        currentGroup = set([current])

        
    
        def recursionSadge(curr, start=1):
            if curr == bito:
                return True
            for i in range(start, n-1):
                add = curr ^^ switchMap[i]
                currentGroup.add(add)
                if recursionSadge(add, i+1 ) == True:
                    return True
            return False
        
        return recursionSadge(current)
        
    def __and__(self, other):
        #Returns true if the two firingclasses have the same supercritical lists. Overrides fc1 & fc2
        if isinstance(other,FiringClass):
            ret = True
            ret = ret and (sorted(self.supers) == sorted(other.supers)) and (sorted(self.supers_signed) == sorted(other.supers_signed)) and (sorted(self.supers_preimage) == sorted(other.supers_preimage))
            ret = ret and (sorted(self.criticals) == sorted(other.criticals)) and (sorted(self.criticals_signed) == sorted(other.criticals_signed)) and (sorted(self.criticals_preimage) == sorted(other.criticals_preimage))
            return ret
            
        return False
    
    def __eq__(self, other):
        #Returns true if fc1 and fc2 have the same L,G. Overrides fc1 == fc2
        if isinstance(other,FiringClass):
            return self.L == other.L and self.G == other.G
        return False
        
    def __ne__(self, other):
        #Returns true if fc1 and fc2 do not have the same L,G. Overrides fc1 != fc2
        return not self.__eq__(other)


    def ones(self):
        return vector([1]*(self.n-1))

    def zeros(self):
        return vector([0]*(self.n-1))

    def getSinkAdjacents(self):
        #Returns the binary vector with 1 on sink-adjacent vertices.
        return vector(sum(self.M))
        
            
        
        