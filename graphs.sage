#GRAPH CREATION TOOLS
def flipEdges(L, lists):
    '''
    Flips specified edges of the graph given by L and returns the new laplacian.
    
    Arguments:
    L     -- A Laplacian matrix to be switched
    lists -- The list of edges to be made negative in a laplacian matrix. Use 1 based indexing!
    '''
    L=copy(L)
    for a,b in lists:
        L[a-1,b-1] *=-1
        L[b-1,a-1] *=-1
    return L

def flipAllEdges(L):
    '''
    Flips all edges of the graph given by L and returns the new laplacian.
    
    Arguments:
    L -- A Laplacian matrix to be switched
    '''
    L = copy(L)
    for a,b in Combinations(L.ncols(),2):
        L[a-1,b-1] *= -1
        L[b-1,a-1] *= -1
    return L

def flipRandomEdges(L, p=0.5):
    '''
    Flips random edges of the graph given by L and returns the new laplacian.
    
    Arguments:
    L -- A Laplacian matrix
    p -- The probability [in RR, 0-1] that any edge is flipped
    '''
    L = copy(L)
    for a,b in Combinations(L.nrows(),2):
        if(p>=random.uniform(0,1)):
            L[a-1,b-1] *= -1
            L[b-1,a-1] *= -1
    return L

def vertexSwitch(L, v):
    '''
    Flips all the edges incident to v in the graph given by L
    
    Arguments:
    L -- A Laplacian matrix
    v -- An integer representing which vertex to switch
    '''
    L = copy(L)
    
    if v>L.nrows():
        return L
    
    for i in range(L.nrows()):
        L[v-1,i]*=-1
        L[i,v-1]*=-1
    return L

def vertexSwitches(L,V):
    '''
    Uses vertexSwitch(L,v) for all the vertices in V
    
    Arguments:
    L -- A Laplacian matrix
    V -- A list of integers representing vertices to switch
    '''
    
    L = copy(L)
    for v in V:
        L = vertexSwitch(L,v)
    return L

def vertexSwitchRandom(L):
    '''
    Switches a sequence of random vertices in the graph given by L
    
    Arguments:
    L -- A Laplacian matrix
    '''
    
    switches = [random.randint(1,L.ncols()) for i in range(L.ncols()-1)]
    return vertexSwitches(L, switches)

def projectFlips(L,M):
    '''
    Flips the edges of M that are flipped in L

    Arguments:
    L -- A laplacian matrix
    M -- Another laplacian matrix of the same size as L
    '''
    L = copy(L)
    M = copy(M)
    for i in range(L.nrows()):
        for j in range(i+1, L.nrows()):
            if M[i,j] != 0 and L[i,j] != 0:
                M[i,j] = L[i,j]
    return M

#GRAPH CREATION TOOLS

def MakeGraph(edges):
    '''
    Creates a graph object and its corresponding laplacian by specifying the edges.
    
    Arguments:
    edges -- A list of edges to include in the graph.
    '''
    G = Graph(edges)
    L = G.laplacian_matrix()[:-1,:-1]
    return L, G

def MakeFromLaplacian(L):
    '''
    Creates a graph object from a reduced laplacian matrix.

    Arguments:
    L -- A reduced laplacian matrix.
    '''
    L = copy(L)
    degtemp = 0
    edges = []
    for i in range(L.ncols()):
        degtemp = 0
        for j in range(L.nrows()):
            if L[i,j] != 0 and i!= j:
                if j>i:
                    edges += [[i+1,j+1]]
                degtemp += 1
        if L[i,i] != degtemp:
            edges += [[i+1,L.ncols()+1]]
                
    G = Graph(edges)
    return L, G


def SubGraph(L, G, vertices):
    '''
    Creates the subgraph of G using only the vertices specified in vertices
    
    Arguments:
    L -- The Laplacian matrix of G
    G -- A graph object
    vertices -- A list of vertices to include in the subgraph G'
    '''
    vertices = copy(vertices)
    L = copy(L)
    G = copy(G)
    
    if not len(G) in vertices:
        vertices += [len(G)]
    
    vertices.sort()
    vertices = vertices[:vertices.index(len(G))+1]
    G = G.subgraph(vertices)
    G.relabel({vertices[i]:i+1 for i in range(len(vertices))})
    if not G.is_connected():
        raise Exception("Subgraph is disconnected")
    L = L.delete_rows([i for i in range(L.ncols()) if not i+1 in vertices]).delete_columns([i for i in range(L.ncols()) if not i+1 in vertices])
    
    Li = G.laplacian_matrix()[:-1,:-1]
    for i in range(Li.ncols()):
        for j in range(i+1, Li.ncols()):
            Li[i,j] = L[i,j]
            Li[j,i] = L[j,i]
    return Li, G

def SubGraphRemoveEdges(L,G, toRemove):
    L = copy(L)
    G = copy(G)
    for edge in toRemove:
        L[edge[0]-1, edge[1]-1] = 0
        L[edge[1]-1, edge[0]-1] = 0
    return MakeFromLaplacian(L)
    
def RandomGraph(n=-1,e=2):
    '''
    Creates a random graph with specified conditions.
    
    Arguments:
    n -- Number of vertices on G. By default gives a random number between 5 and 9
    e -- The average number of edges. By default gives 2. Note that passing e=1 generates a random tree.
    '''
    if n == -1:
        n = randint(5,9)
        
    if e == 1:
        G = graphs.RandomBarabasiAlbert(n,e)
    else:
        e = 0.5*(e+1)*n
        G = graphs.RandomGNM(n, e)
        while not G.is_connected():
            G = graphs.RandomGNM(n, e)
    
    G.relabel({i:i+1 for i in range(n)})
    
    L = G.laplacian_matrix()[:-1,:-1]
    return L, G

def RandomSpecGraph(n,e):
    '''
    Creates a random connected graph with n nodes and e edges.
    Throws an error if e<n-1.

    Arguments:
    n -- Number of vertices on G
    e -- Number of edges on G
    '''
    if e<n-1:
        raise Exception("Disconnected Graph, not enough edges in RandomSpecGraph()")

    G = graphs.RandomBarabasiAlbert(n,1)
    G.relabel({i:i+1 for i in range(n)})
    e -= n-1
    while e>0:
        v1 = randint(1,n)
        v2 = randint(1,n)
        if v1 != v2 and not G.has_edge([v1,v2]):
            G.add_edge([v1,v2])
            e-=1

    return G.laplacian_matrix()[:-1,:-1], G
    
def RandomPGraph(n,p, connected = True):
    '''
    Creates a random graph with n nodes and probability p of an edge between v1 and v2

    Arguments:
    n -- The number of vertices on G
    p -- The probability of an edge between v1 and v2
    '''
    if connected and p< (2/n):
        raise Exception("Cannot make graph connected")

    G = randomGNP(n,p)

    if connected:
        while not G.connected():
            G = randomGNP(n,p)
    return G

# Classes of graphs

def Cycle(n=8):
    '''
    Creates a graph object representing the cyclic graph Cn
    
    Arguments:
    n -- The number of vertices in Cn including the sink. 
    '''
    G = Graph([[a,b] for a,b in pairwise(range(1, n+1))] + [[1,n]])
    L = G.laplacian_matrix()[:-1,:-1]
    return L, G

def Complete(n=6):
    '''
    Creates a graph object representing the complete graph Kn
    
    Arguments:
    n -- The number of vertices in Kn including the sink
    '''
    G = Graph([[a,b] for a,b in combinations(range(1,n+1),2)])
    L = G.laplacian_matrix()[:-1,:-1]
    return L, G

def Cube(n=3):
    '''
    Creates a graph object representing the cube graph Qn
    
    Arguments:
    n -- The dimension of the cube. Note that Qn will have 2^n vertices.
    '''
    G = graphs.CubeGraph(n)
    G.relabel({a:i+1 for a,i in zip(sorted(list(G.get_vertices().keys())), range(2^n))})
    
    L = G.laplacian_matrix()[:-1, :-1]
    return L,G

def Bipartite(m=3,n=4):
    '''
    Creates a graph object representing the complete bipartite graph Km,n
    
    Arguments:
    m -- The number of vertices in one partition
    n -- The number of vertices in the other partition
    '''
    G = Graph([[a+1,b+m+1] for a,b in product(range(m),range(n))])
    L = G.laplacian_matrix()[:-1,:-1]
    return L,G

def Wheel(n=5):
    '''
    Creates a graph object representing the wheel graph Wn. Note that the sink will always be the middle vertex.
    
    Arguments:
    n -- The number of vertices on the exterior cycle. Note that Wn has n+1 vertices
    '''
    cycle = [[a,b] for a,b in pairwise(range(1,n+1))] + [[1,n]]
    spokes = [[a, n+1] for a in range(1, n+1)]
    G = Graph(cycle + spokes)
    L = G.laplacian_matrix()[:-1,:-1]
    return L, G
    
def Fan(n=6):
    '''
    Creates a graph object representing the fan graph Fn. Note that the sink will always be the fanning vertex.
    
    Arguments:
    n -- The number of vertices in Fn.
    '''
    cycle = [[a,b] for a,b in pairwise(range(1,n+1))]
    spokes = [[a, n] for a in range(1, n - 1)]
    G = Graph(cycle + spokes)
    L = G.laplacian_matrix()[:-1,:-1]
    return L, G

def FavoriteGraph(wrong = False):
    '''
    Creates a graph object representing Kyle, the complete graph minus a non-sink edge.
    '''
    if wrong:
        G = Graph([[1,2],[2,3],[4,3],[4,1],[2,4]])
    else:
        G = Graph([[1,2],[2,3],[4,3],[4,1],[1,3]])
    L = G.laplacian_matrix()[:-1,:-1]
    return L, G

def displayGraph(L,G):
    '''
    Displays the signed graph given by L, G visually.
    
    Arguments:
    L -- The Laplacian matrix of G
    G -- A graph object
    '''
    G = copy(G)
    L = copy(L)
    n = L.nrows()
    graphstr = "M, G = MakeGraph(["
    edgestr = "L = flipEdges(M, ["
    edgestr2 = ""
    for a ,b , c in G.edges():
        G.set_edge_label(a,b,'+')
        graphstr += f"[{a},{b}], "
        if a-1 != n and b-1 != n and L[a-1,b-1] == 1:
            G.set_edge_label(a,b,'–')
            edgestr2 += f"[{a},{b}], "
    
    graphstr = graphstr[:-2] + "])"
    edgestr += edgestr2[:-2] + "])\n"
    if edgestr2 == "":
        edgestr = ""
    
    G.relabel({L.ncols()+1: 'q'})
    G.show(edge_labels=True, vertex_colors = {'magenta':['q']}, color_by_label={'–':'red'})
    print(f"Generated by: \n{graphstr}\n{edgestr}")
    
    
    