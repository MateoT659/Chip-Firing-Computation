#HELPER FUNCTIONS

def timeSince(timeT, startFrom = False, label = "Time"):
    '''
    Prints the time elapsed since timeT was defined in HH, MM, SS.SSS format
    
    Arguments:
    timeT     -- The time given by time.time() established prior to this function's calling
    startFrom -- If True, return the current time to be stored in a new timeT variable.
    label     -- The label that will be printed with the time. 
    '''
    t = time.time() - timeT
    s = t%60
    t = int(t)
    m = t//60
    h = m//60
    print(f"{label}: {int(h):02d} h, {int(m):02d} m, {s:02.3f} s")
    if startFrom:
        return time.time()

def displaySwitchClasses(G, full = False, show = True, metric = 'positive'):
    '''
    Displays the structure of the switching equivalence classes of G along with representatives from each class.
    
    Arguments:
    G -- A graph object
    full -- When True, displays a graph from every equivalence class. Otherwise, only one graph is shown for each unique determinant.
    show -- When True, displays the graphs visually. Otherwise, only display the group structure and determinants.
    metric -- The sorting metric of the graph representative. 'positive'/'negative' returns the graph with the most/least positive edges
              'posVerts' returns the graph with the most positive incident vertices
    '''
    classes = getSwitchClasses(G,metric)
    fc = FiringClass(classes[0],G)
    fc.displayShallow()

    firings = []
    for Li in classes:
        fc = FiringClass(Li,G)
        firings += [fc]
    firings.sort(key = lambda x: x.dL)

    dLmap = {}
    for f in firings:
        dLmap[f.dL] = dLmap[f.dL]+1 if f.dL in dLmap else 1

    print(f"\nNumber of different determinants: {len(dLmap)}")
    dL = 0
    determinants = []
    for a in firings:
        if dL != a.dL:
            dL = a.dL
            determinants += [dL]
            print(f"\n\n Determinant {a.dL}, {dLmap[a.dL]} class{'es' if dLmap[a.dL] != 1 else ''} \n")
            if not full:
                if show:
                    a.displayGraph()
                a.displayGroupStructure()
        if full:
            if show:
                a.displayGraph()
            a.displayGroupStructure()
            
    return determinants
    
    
def columnDisplay(lists, labels = [], limit = False):
    '''
    Displays lists next to each other in a column format.
    
    Arguments:
    lists  -- A list of lists to display in order
    labels -- A list of labels for each list. Defaults to "List i"
    limit  -- When True, suppresses the output based on MAX_CONF_DISPLAY
    '''
    maxLen = max([len(A) for A in lists])
    
    if limit and maxLen>MAX_CONF_DISPLAY:
        print("[Suppressed - Too many entries!]")
        return
    
    lists = [A+[None]*(maxLen - len(A)) for A in lists]
    
    maxStrLens = [max([len(str(a)) for a in A]) for A in lists]
    
    for i in range(len(labels), len(lists)):
        labels.append(f"List {i+1}")
    
    zipped = [tuple([a[i] for a in lists]) for i in range(maxLen)]

    for i in range(maxLen):
        pstr = ""
        for j in range(len(labels)):
            pstr += f"{labels[j]}: {zipped[i][j]}{' '*(maxStrLens[j] - len(str(zipped[i][j])))}\t"
        print(pstr)
    
    
def disjoint(L1, L2):
    #Returns the disjoint of lists L1 and L2
    return [copy(l) for l in L1 if l in L2]


def disjointAll(L):
    #Returns the disjoint of all lists in L
    leng = len(L)
    if leng == 0:
        return []
    if leng == 1:
        return L[0]
    if leng == 2:
        return disjoint(L[0], L[1])
    else:
        return disjoint(disjointAll(L[:leng//2]), disjointAll(L[leng//2:]))

def disjointUnion(L):
    '''
    Returns any elements that are in more than one of the lists in L
    DU = (L1 and L2) or (L1 and L3) or ... (Ln-1 and Ln)
    Arguments:
    L -- A list of lists
    '''
    ret = []
    for i in range(len(L)):
        for j in range(i+1, len(L)):
            ret += [disjoint(L[i], L[j])]
    return ret
        
def get_sp_config(sp,v):
    '''
    Returns the sandpileConfig object of a vector v
    
    Arguments:
    sp -- A sandpile object
    v  -- A vector or list
    '''
    return SandpileConfig(sp, list(v)) 

def equivalentSuperstable(sp, vec):
    return vector(get_sp_config(sp,vec).equivalent_superstable(False).values())

def equivalentCritical(sp,vec):
    return vector(get_sp_config(sp,vec).equivalent_recurrent(False).values())

def list_subtract(L1, L2):
    #Returns the list given by L1 - L2
    return [l for l in L1 if l not in L2]

def floor_vec(v):
    #Returns floor(v) for a vector v
    v = copy(v)
    return vector([floor(i) for i in v])

def ciel_vec(v):
    #Returns cieling(v) for a vector v
    v = copy(v)
    return vector([ceil(i) for i in v])


def frac_vec(v):
    #Returns v - floor(v) for a vector v
    v = copy(v)
    return v - floor_vec(v)

def modQQ(a,m):
    b = a%m
    return (b + m) if b<0 else b

def mod_vec(v, m):
    #Returns v mod m for a vector v and real m
    v = copy(v)
    return vector([modQQ(i,m) for i in v])

def abs_vec(v):
    #Returns absolute_value(v) for a vector v
    v = copy(v)
    return vector([abs(i) for i in v])


def abs_matrix(A):
    #Returns absolute_value(A) for a matrix A
    A = copy(A)
    for i in range(A.ncols()):
        for j in range(A.nrows()):
            A[j,i] = abs(A[j,i])
    return A


def integer_list(L):
    #Checks if a list is composed of only integers
    for a in L:
        if not (a in ZZ):
            return False
    return True


def get_cmax(L):
    #Returns maximal critical vector of a Laplcian matrix L
    return  vector([L[a][a]-1 for a in range(L.nrows())])


def sstab(c):
    #Returns the equivalent superstable of the rational vector c
    fl = SandpileConfig(c.sandpile(), {k: floor(v) for k,v in c.items()})
    return fl.equivalent_superstable() + c - fl


def find_order(c):
    #Finds the order of c in graph G
    copy = sstab(c+c)
    order = 1
    while copy != sstab(c):
        copy = sstab(copy+c)
        order+=1
    return order


def cycleSum(v):
    #Returns the cycle sum used in the proofs of cyclic graphs
    return sum([ (-1)^(i+2) * (i+1) * x for i,x in enumerate(v)])%(len(v)+1)

def value_counts(list):
    '''
    Returns a value count dictionary of the list given, with keys being elements of list. 
    Elements must be hashable; use value_counts vec for vectors

    Arguments:
    list -- A list of hashables to be counted
    '''
    d = {}
    for v in list:
        d[v] = d[v]+1 if v in d else 1
    return d

def value_counts_vec(vectors):
    '''
    Returns a value count dictionary of the vectors in vectors using tuples as keys.

    Arguments:
    vectors -- A list of vectors or lists to be counted.
    '''

    tuples = [tuple(v) for v in vectors]
    return value_counts(tuples)

def count_sublists(list_of_lists):
    #Returns a dictionary containing the amount of times a sublist is in a list of lists, indexed by sublist
    sublist_count = {}
    for sublist in list_of_lists:
        sublist_tuple = tuple(sublist)
        if sublist_tuple in sublist_count:
            sublist_count[sublist_tuple] += 1
        else:
            sublist_count[sublist_tuple] = 1
    return sublist_count


def primaryDecomposition(L):
    #Finds the primary decomposition of a signed laplacian L
    L=copy(L)
    ed = L.elementary_divisors()
    A = [ed[i] for i in range(L.nrows()) if ed[i]>1]
    return len(set(A))

def invariantFactors(L):
    #Finds the invariant factors of a signed laplacian L
    L = copy(L)
    ed = L.elementary_divisors()
    return [ed[i] for i in range(L.nrows()) if ed[i]>1]



def getGroupStructure(L):
    '''
    Calculates the group structure of the graph given by L and gives a string showing it as the direct product of integers
    
    Arguments:
    L - A signed Laplacian matrix
    '''
    ed = L.elementary_divisors()
    A = [i for i in ed if i > 1]
    ret = ""

    for a in A:
        ret += f"Z{a} ⊕ "
    ret = ret[:-3]
    return ret


def numSwitchClasses(G):
    #Gives the number of switching classes in a graph G. if G is disconnected or E<V, may return rationals.
    G = copy(G)
    G.delete_vertex(G.num_verts())
    return 2^(G.num_edges() - G.num_verts() + 1)


def sizeSwitchClasses(G):
    #Gives the size of each switching class in G. Each class has the same size.
    G = copy(G)
    G.delete_vertex(G.num_verts())
    return 2^(G.num_verts()-1)


def numPositiveVerts(L, G):
    '''
    Finds the number of vertices with all positive incident edges in the signed graph L,G
    
    Arguments:
    L -- The Laplacian matrix of G
    G -- A graph object
    '''
    L = copy(L)
    G = copy(G)
    M = G.laplacian_matrix()[:-1,:-1]
            
    LMinv = L/M            
    I = matrix.identity(L.nrows())
    count = 0
    for i in range(LMinv.nrows()):
        if LMinv[i] == I[i]:
            count+=1
    return count


def numNegativeEdges(L):
    '''
    Finds the number of negative edges in the signed graph given by L
    
    Arguments:
    L -- A Laplacian matrix
    '''
    L = copy(L)
    c = 0
    for i in range(L.nrows()):
        for j in range(i+1,L.ncols()):
            if L[i,j] == 1:
                c += 1
    return c

def getSwitchClasses(G, metric = 'positive'):
    '''
    Returns a list of laplacian matrices representing signed graphs from each switching equivalence class of G
    
    Arguments:
    G      -- A graph object
    metric -- The metric through which the representatives will be picked.
              'positive'/'negative' finds the Laplacian with the most/least positive edges
              'numVerts' finds the Laplacian with the most positive incident vertices
    '''
    
    #returns list of laplacian graphs representing each equivalence class. 
    G = copy(G)
    L = G.laplacian_matrix()[:-1,:-1]
    #randomind2^n == random bitset
    ret = []
    seen = set()
    n = G.num_verts()
    nClasses = numSwitchClasses(G)
    sizeClasses = sizeSwitchClasses(G)
    nStates = nClasses*sizeClasses
    #indexMap maps which index of a bitstring corresponds to which edge
    indexMap = [] #indexes of MATRIX
    switchMap = [Bitset() for i in range(n)] #switchMap[0] = 0
    for i in range(n-1):
        for j in range(i+1, n-1):
            if L[i,j] != 0:
                indexMap += [(i,j)]
                switchMap[i+1].add(len(indexMap)-1)
                switchMap[j+1].add(len(indexMap)-1)     
    
    l = len(indexMap)
        
    current = FrozenBitset(capacity=l)
    
    def recursionSadge(curr, start=1):
        for i in range(start, n-1):
            add = curr ^^ switchMap[i]
            currentGroup.add(add)
            recursionSadge(add, i+1)
    x = 0
    while len(seen) < nStates:
        while current in seen:
            #try random bitsets until you find new one
        
            current = FrozenBitset([i for i in range(l) if x & (2^i) != 0], capacity = l)
            x += 1
            
        currentGroup = set([current])
        recursionSadge(current)
        
        seen |= currentGroup
        if metric == 'negative':
            ret += [max(list(currentGroup), key = lambda x:(len(x),x))]
        elif metric == 'positive':
            ret += [max(list(currentGroup), key = lambda x:(-len(x),x))]
        else:
            ret += [max(list(currentGroup), key = lambda x: sum([1 for i in range(1,n) if (x&switchMap[i]) == switchMap[0] ]))]
    retM = []
        
    for bits in ret:
        Li = copy(L)
        for bit in bits:
            Li[indexMap[bit]] = 1
            Li[indexMap[bit][::-1]] = 1
        retM += [Li]
    return retM


def getSwitchClass(L,G,metric='posVerts'):
    '''
    Returns a list of all Laplacian matrices switching equivalent L, sorted by the metric
    '''
    G = copy(G)
    L = copy(L)
    sizeClass = sizeSwitchClasses(G)
    switchClass = set()
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
    
    currentGroup = set([current])
    
    #extremely mid solution
#     while len(currentGroup) < sizeClass:
#         toAdd = set()
#         for bits in currentGroup:
#             for i in range(1,n):
#                 add = bits^^switchMap[i]
#                 if add not in currentGroup:
#                     toAdd.add(bits^^switchMap[i])
#                     pathdict[add] = tuple(list(pathdict[bits]) + [i])
#         currentGroup |= toAdd
    validinds = [i+1 for i in range(n-2)]
    
    #brute force
#     for mask in product([0,1], repeat = n-2):
#         add = current
#         path = []
#         for i in compress(validinds, mask):
#             add^^= switchMap[i]
#             path+=[i]
#         currentGroup.add(add)
    
    #recursive, uses only one xor for each guy. next goal make this not recursive
    def recursionSadge(curr, start=1):
        for i in range(start, n-1):
            add = curr ^^ switchMap[i]
            currentGroup.add(add)
            recursionSadge(add, i+1)
    
    recursionSadge(current)
        
    if metric == 'negative':
        ret = sorted(list(currentGroup), key = lambda x:(len(x),x))
    elif metric == 'positive':
        ret = sorted(list(currentGroup), key = lambda x:(-len(x),x))
    else:
        ret = sorted(list(currentGroup), key = lambda x: sum([1 for i in range(1,n) if (x&switchMap[i]) == switchMap[0] ]))

    lapret = []
    for bits in ret:
        Li = G.laplacian_matrix()[:-1,:-1]
        for bit in bits:
            Li[indexMap[bit]] = 1
            Li[indexMap[bit][::-1]] = 1
        lapret += [Li]
        
    return lapret

def getSwitchRep(L, G, metric='posVerts'):
    '''
    Returns the Laplacian matrix switching equivalent to L that best fits the metric
    
    Arguments:
    L -- The Laplacian matrix given by G
    G      -- A graph object
    metric -- The metric through which the representatives will be picked.
              'positive'/'negative' finds the Laplacian with the most/least positive edges
              'numVerts' finds the Laplacian with the most positive incident vertices
    '''
    
    G = copy(G)
    L = copy(L)
    sizeClass = sizeSwitchClasses(G)
    switchClass = set()
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
    
    currentGroup = set([current])
    
    #extremely mid solution
#     while len(currentGroup) < sizeClass:
#         toAdd = set()
#         for bits in currentGroup:
#             for i in range(1,n):
#                 add = bits^^switchMap[i]
#                 if add not in currentGroup:
#                     toAdd.add(bits^^switchMap[i])
#                     pathdict[add] = tuple(list(pathdict[bits]) + [i])
#         currentGroup |= toAdd
    validinds = [i+1 for i in range(n-2)]
    
    #brute force
#     for mask in product([0,1], repeat = n-2):
#         add = current
#         path = []
#         for i in compress(validinds, mask):
#             add^^= switchMap[i]
#             path+=[i]
#         currentGroup.add(add)
    
    #recursive, uses only one xor for each guy. next goal make this not recursive
    def recursionSadge(curr, start=1, past=[]):
        for i in range(start, n-1):
            add = curr ^^ switchMap[i]
            currentGroup.add(add)
            recursionSadge(add, i+1,past+[i] )
    
    recursionSadge(current)
        
    if metric == 'negative':
        ret = max(list(currentGroup), key = lambda x:(len(x),x))
    elif metric == 'positive':
        ret = max(list(currentGroup), key = lambda x:(-len(x),x))
    else:
        ret = max(list(currentGroup), key = lambda x: sum([1 for i in range(1,n) if (x&switchMap[i]) == switchMap[0] ]))
        
    Li = G.laplacian_matrix()[:-1,:-1]
    for bit in ret:
        Li[indexMap[bit]] = 1
        Li[indexMap[bit][::-1]] = 1
        
    numMap = {}
    path = getSwitches(L,Li,G)
    
    for p in path:
        numMap[p] = numMap[p]+1 if p in numMap else 1
    path = []
    
    for k in numMap:
        if numMap[k]%2 == 1:
            path += [k]
        
    return Li, path

def switchingEquivalent(L,M):
    '''
    Determines if the two matrices L and M are switching equivalent

    Arguments:
    L -- A Laplacian matrix of G
    M -- Another Laplacian matrix of G
    '''
    L = copy(L)
    M = copy(M)
    n = L.nrows()+1
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
        if(M[indexMap[i]] == 1):
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
    

def getSwitches(L,M,G):
    '''
    Finds the sequence of vertices that must be switched to convert between L and M
    
    Arguments:
    L -- A Laplacian matrix of G
    M -- Another Laplacian matrix of G that is switching equivalent to L
    G -- A graph object
    '''
    G = copy(G)
    L = copy(L)
    M = copy(M)
    swap = L-M
    queue = []
    seen = []
    path = []
    n = M.nrows()
    G.delete_vertex(n+1)
    Gs = G.connected_components_subgraphs()
    for G in Gs:
        queue = [G.random_vertex()]
        while queue:
            curr = queue.pop(-1)
            for i in range(n):
                if swap[curr-1, i] != 0:
                    L = vertexSwitch(L, i+1)
                    swap = L-M
                    path += [i+1]
        
            seen += [curr]
        
            queue = list_subtract(G.neighbors(curr), seen+queue) + queue
    
    return tuple(path)

def getCycleSwitches(L,M):
    '''
    Finds the sequence of vertices that must be switched to convert between L and M for cyclic graphs.
    
    Arguments:
    L -- A Laplacian matrix of a cyclic graph
    M -- Another Laplacian matrix switching equivalent to L
    '''
    
    #scan for 4s, then scan for 2s
    L = copy(L)
    M = copy(M)
    swap = (L-M)/2
    path = []
    n = L.nrows()
    zeros = matrix(n, n)
    while swap != zeros:
        #check for 4s and swap them (not outsides)
        for i in range(1, n-1):
            if swap[i,i-1] == 1 and swap[i,i+1] == 1 and swap[i-1,i] == 1 and swap[i+1, i] == 1:
                swap[i,i-1] = swap[i,i+1] = swap[i+1,i] = swap[i-1, i] = 0
                path += [i+1]
        #check for 2s and swap them
        for i in range(0,floor(n/2)):
            if swap[i,i+1] == 1 and swap[i+1, i] == 1:
                swap[i,i+1] = swap[i+1,i] = 0
                if i != 0:
                    swap[i-1,i] = swap[i,i-1] = 1
                path+= [i+1]
        for i in range(floor(n/2), n):
            if swap[i,i-1] == 1 and swap[i-1,i] == 1:
                swap[i,i-1] = swap[i-1,i] = 0
                if i != n-1:
                    swap[i,i+1] = swap[i+1,i] = 1
                path += [i+1]
                
    # only keep odds
    numMap = {}
    for p in path:
        numMap[p] = numMap[p]+1 if p in numMap else 1
    path = []
    for k in numMap:
        if numMap[k]%2 == 1:
            path += [k]
    return tuple(path)


def isCycle(M):
    #Checks if the unsigned Laplacian M represents a cyclic graph 
    n = M.nrows()
    L, G = Cycle(n+1)
    
    return L-M == matrix(n,n)


def tupleM(M):
    '''
    Converts a matrix M into a hashable tuple matrix. Convert back using M = matrix(tupleM(M))
    
    Arguments:
    M - A matrix
    '''
    return tuple([tuple(m) for m in M])

def dhar_burn(G, c):
    '''
    Performs Dhar's Burning Algorithm on the graph G and configuration c
    Determines whether the configuration c is superstable or not.

    Arguments:
    G -- A graph object
    c -- A chip configuration
    '''
    d = len(c)
    burning = [0] * d
    goal = [-1]*d
    queue = [d+1]
    curr = -1
    while len(queue)>0:
        curr = queue.pop(-1)
        for n in G.neighbors(curr):
            if n<=d and burning[n-1] >=0:
                burning[n-1] += 1
                if burning[n-1]>c[n-1]:
                    queue = [n]+queue
                    burning[n-1] = -1
                            
    return vector(burning) == vector(goal)

def le_vec(v1, v2):
    '''
    Checks if the vector v1 is elementwise greater than v2

    Arguments
    v1 -- A vector of length d
    v2 -- Another vector of length d
    '''
    if len(v1) != len(v2):
        return False
    
    for i in range(len(v1)):
        if v2[i]<v1[i]:
            return False
        
    return True

def ge_vec(v1, v2):
    '''
    Checks if the vector v1 is elementwise greater than v2

    Arguments
    v1 -- A vector of length d
    v2 -- Another vector of length d
    '''
    if len(v1) != len(v2):
        return False
    
    for i in range(len(v1)):
        if v2[i]>v1[i]:
            return False
        
    return True

def positive_vec(v):
    #Checks if the vertex v has all non-negative entries
    for a in v:
        if a<0:
            return False
    return True

def list_prod(L):
    p = 1
    for l in L:
        p*=l
    return p

def first(iterable):
    #Return the first element of an iterable (list, set, dict, etc)
    for e in iterable:
        break
    return e


def check_dominant(L1,L2):
    #Check if a list has larger entries than another
    return not False in [x <= y for x,y in zip(L1,L2)]


def collect_minimals(L):
    return [v for v in L if not True in [check_dominant(w,v) for w in L if not w ==v]]
