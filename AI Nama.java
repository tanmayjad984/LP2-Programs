LP 2 


AI Problem Statements

######  1) Implement depth first search algorithm and Breadth First Search
algorithm, use an undirected graph and develop a recursive algorithm for
searching all the vertices of a graph or tree data structure. ######


#Graph.java


import java.util.*;

public class Graph {
    private int V; // number of vertices
    private LinkedList<Integer>[] adj; // adjacency list

    // constructor
    public Graph(int V) {
        this.V = V;
        adj = new LinkedList[V];
        for (int i = 0; i < V; i++) {
            adj[i] = new LinkedList<>();
        }

    }

    // add edge to the graph
    public void addEdge(int v, int w) {
        adj[v].add(w);
        adj[w].add(v);
    }

    // Depth-first search
    public void dfs(int v) {
        boolean[] visited = new boolean[V];
        dfsUtil(v, visited);
    }

    private void dfsUtil(int v, boolean[] visited) {
        visited[v] = true;
        System.out.print(v + " ");

        Iterator<Integer> it = adj[v].listIterator();
        while (it.hasNext()) {
            int n = it.next();
            if (!visited[n]) {
                dfsUtil(n, visited);
            }
        }
    }

    // Breadth-first search
    public void bfs(int v) {
        boolean[] visited = new boolean[V];
        LinkedList<Integer> queue = new LinkedList<>();
        visited[v] = true;
        queue.add(v);

        while (!queue.isEmpty()) {
            int n = queue.poll();
            System.out.print(n + " ");

            Iterator<Integer> it = adj[n].listIterator();
            while (it.hasNext()) {
                int m = it.next();
                if (!visited[m]) {
                    visited[m] = true;
                    queue.add(m);
                }
            }
        }
    }

    public static void main(String[] args) {
        Graph g = new Graph(6);
        g.addEdge(0, 1);
        g.addEdge(0, 2);
        g.addEdge(1, 3);
        g.addEdge(2, 4);
        g.addEdge(3, 4);
        g.addEdge(3, 5);
    
        System.out.print("DFS: ");
        g.dfs(0);
        System.out.println();
    
        System.out.print("BFS: ");
        g.bfs(0);
        System.out.println();
    
    
    }
}



######2) Implement A Star Algorithm for 8 puzzle problem######

#AStarPuzzle.java

import java.util.*;

public class AStarPuzzle {

    private static final int[] dx = {-1, 0, 1, 0};
    private static final int[] dy = {0, -1, 0, 1};

    public static void main(String[] args) {
        int[][] initialState = {{1, 2, 3}, {4, 0, 5}, {6, 7, 8}};
        int[][] goalState = {{1, 2, 3}, {4, 5, 6}, {7, 8, 0}};

        // Create a priority queue to store the open nodes.
        PriorityQueue<Node> open = new PriorityQueue<>((a, b) -> (a.cost + a.h) - (b.cost + b.h));

        // Add the initial state to the open list.
        open.add(new Node(initialState, 0, calculateHScore(initialState, goalState), null));

        while (!open.isEmpty()) {
            // Get the node with the lowest f-score.
            Node current = open.poll();

            // If the current node is the goal state, then we have found a solution.
            if (Arrays.deepEquals(current.state, goalState)) {
                System.out.println("Found a solution!");
                printPath(current);
                return;
            }

            // For each of the neighbors of the current node,
            for (int i = 0; i < 4; i++) {
                int x = current.getBlankX() + dx[i];
                int y = current.getBlankY() + dy[i];

                // If the neighbor is within the bounds of the board,
                if (x >= 0 && x < 3 && y >= 0 && y < 3) {
                    // And the neighbor is not the empty space,
                    int value = current.state[x][y];
                    if (value != 0) {
                        // Then create a new node with the neighbor's state.
                        int[][] newState = deepCopy(current.state);
                        newState[x][y] = 0;
                        newState[current.getBlankX()][current.getBlankY()] = value;

                        // Calculate the h-score for the neighbor.
                        int h = calculateHScore(newState, goalState);

                        // Create the neighbor node and add it to the open list.
                        Node neighbor = new Node(newState, current.cost + 1, h, current);
                        open.add(neighbor);
                    }
                }
            }
        }

        // If we reach this point, then there is no solution.
        System.out.println("No solution found!");
    }

    private static int calculateHScore(int[][] state, int[][] goalState) {
        int h = 0;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (state[i][j] != goalState[i][j]) {
                    h++;
                }
            }
        }
        return h;
    }

    private static void printPath(Node node) {
        List<Node> path = new ArrayList<>();
        while (node != null) {
            path.add(node);
            node = node.parent;
        }

        Collections.reverse(path);

        for (Node n : path) {
            System.out.println(Arrays.deepToString(n.state));
        }
    }

    private static int[][] deepCopy(int[][] array) {
        int[][] copy = new int[3][3];
        for (int i = 0; i < 3; i++) {
            System.arraycopy(array[i], 0, copy[i], 0, 3);
        }
        return copy;
    }
}

class Node {
    int[][] state;
    int cost;
    int h;
    Node parent;

    public Node(int[][] state, int cost, int h, Node parent) {
        this.state = state;
        this.cost = cost;
        this.h = h;
        this.parent = parent;
    }

    public int getBlankX() {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (state[i][j] == 0) {
                    return i;
                }
            }
        }
        return -1; // Blank not found
    }

    public int getBlankY() {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (state[i][j] == 0) {
                    return j;
                }
            }
        }
        return -1; // Blank not found
    }
}





###### 3) Implement A Star Algorithm to find out shortest distance ######


#Node.java
import java.util.*;

public class Node {
    private int x;
    private int y;
    private Map<Node, Integer> neighbors; // Neighboring nodes and cost to move to each neighbor
    
    public Node(int x, int y) {
        this.x = x;
        this.y = y;
        this.neighbors = new HashMap<>();
    }
    
    public int getX() {
        return x;
    }
    
    public int getY() {
        return y;
    }
    
    public void addNeighbor(Node neighbor, int cost) {
        neighbors.put(neighbor, cost);
    }
    
    public Collection<Node> getNeighbors() {
        return neighbors.keySet();
    }
    
    public int getCost(Node neighbor) {
        return neighbors.get(neighbor);
    }
    @Override
    public String toString() {
        return "(" + x + "," + y + ")";
    }
}




#AStar.java


import java.util.*;

public class AStar {
    public static int estimateDistance(Node a, Node b) {
        // Heuristic function to estimate distance between nodes a and b
        // In this example, we use the Euclidean distance between the node coordinates
        int dx = a.getX() - b.getX();
        int dy = a.getY() - b.getY();
        return (int) Math.sqrt(dx*dx + dy*dy);
    }

    public static List<Node> findShortestPath(Node start, Node goal) {
        Map<Node, Integer> gScore = new HashMap<>(); // Cost from start along best known path
        gScore.put(start, 0);
        Map<Node, Integer> fScore = new HashMap<>(); // Estimated total cost from start to goal through y
        fScore.put(start, estimateDistance(start, goal));
        PriorityQueue<Node> openSet = new PriorityQueue<>(Comparator.comparingInt(fScore::get));
        openSet.add(start);
        Map<Node, Node> cameFrom = new HashMap<>(); // Parent node on best known path
        while (!openSet.isEmpty()) {
            Node current = openSet.poll();
            if (current.equals(goal)) {
                // Reconstruct the path from start to goal
                List<Node> path = new ArrayList<>();
                path.add(current);
                while (cameFrom.containsKey(current)) {
                    current = cameFrom.get(current);
                    path.add(current);
                }
                Collections.reverse(path);
                return path;
            }
            for (Node neighbor : current.getNeighbors()) {
                int tentativeGScore = gScore.get(current) + current.getCost(neighbor);
                if (!gScore.containsKey(neighbor) || tentativeGScore < gScore.get(neighbor)) {
                    cameFrom.put(neighbor, current);
                    gScore.put(neighbor, tentativeGScore);
                    fScore.put(neighbor, tentativeGScore + estimateDistance(neighbor, goal));
                    if (!openSet.contains(neighbor)) {
                        openSet.add(neighbor);
                    }
                }
            }
        }
        return null; // No path found
    }
    public static void main(String[] args) {
    // Create example graph
    Node start = new Node(0, 0);
    Node a = new Node(1, 2);
    Node b = new Node(3, 5);
    Node c = new Node(5, 2);
    Node goal = new Node(6, 6);
    start.addNeighbor(a, 5);
    start.addNeighbor(c, 10);
    a.addNeighbor(b, 3);
    b.addNeighbor(goal, 2);
    c.addNeighbor(b, 6);
    c.addNeighbor(goal, 4);
    
    // Find shortest path using A* algorithm
    List<Node> path = AStar.findShortestPath(start, goal);
    if (path != null) {
        System.out.println("Shortest path found: " + path);
    } else {
        System.out.println("No path found from start to goal");
    }
}

}






###### 4) Implement a solution for a Constraint Satisfaction Problem using Branch
and Bound and Backtracking for n-queens problem ######


#NQueens.java

public class NQueens {

    private int n;
    private int[] cols;
    private boolean[] usedCols;
    private boolean[] usedDiag1;
    private boolean[] usedDiag2;
    private int solutionsCount;

    public NQueens(int n) {
        this.n = n;
        this.cols = new int[n];
        this.usedCols = new boolean[n];
        this.usedDiag1 = new boolean[2 * n - 1];
        this.usedDiag2 = new boolean[2 * n - 1];
    }

    public void solve() {
        solutionsCount = 0;
        placeQueen(0);
    }

    private void placeQueen(int row) {
        if (row == n) {
            solutionsCount++;
            printSolution();
            return;
        }

        for (int col = 0; col < n; col++) {
            if (isValidPosition(row, col)) {
                cols[row] = col;
                usedCols[col] = true;
                usedDiag1[row + col] = true;
                usedDiag2[row - col + n - 1] = true;

                placeQueen(row + 1);

                cols[row] = 0;
                usedCols[col] = false;
                usedDiag1[row + col] = false;
                usedDiag2[row - col + n - 1] = false;
            }
        }
    }

    private boolean isValidPosition(int row, int col) {
        return !usedCols[col] && !usedDiag1[row + col] && !usedDiag2[row - col + n - 1];
    }

    private void printSolution() {
        System.out.print("Solution " + solutionsCount + ": ");
        for (int i = 0; i < n; i++) {
            System.out.print("(" + i + "," + cols[i] + ") ");
        }
        System.out.println();
    }
    public static void main(String[] args) {
        int n = 4;
        NQueens solver = new NQueens(n);
        solver.solve();
    }
    
}




###### 5) Implement Greedy search algorithm for Selection Sort ######

#selectionSort.java


import java.util.*;
public class selectionSort
{
public static void selectionSort(int[] arr) {
    int n = arr.length;
    for (int i = 0; i < n - 1; i++) {
        int minIndex = i;
        for (int j = i + 1; j < n; j++) {
            // Greedy algorithm to find minimum element
            if (arr[j] < arr[minIndex]) {
                minIndex = j;
            }
        }
        // Swap minimum element with current element
        int temp = arr[i];
        arr[i] = arr[minIndex];
        arr[minIndex] = temp;
    }
}
public static void main(String[] args) {
    int[] arr = {4, 2, 9, 3, 6};
    System.out.println("Unsorted array: " + Arrays.toString(arr));
    selectionSort(arr);
    System.out.println("Sorted array: " + Arrays.toString(arr));
}

}





###### 6) Implement Greedy search algorithm for Prim's Minimal Spanning Tree Algorithm ######



#PrimMST.java


import java.util.*;

public class PrimMST {
    private static final int INF = Integer.MAX_VALUE;

    public static void main(String[] args) {
        int[][] graph = {
                {INF, 2, INF, 6, INF},
                {2, INF, 3, 8, 5},
                {INF, 3, INF, INF, 7},
                {6, 8, INF, INF, 9},
                {INF, 5, 7, 9, INF}
        };

        int[] parent = primMST(graph);

        System.out.println("Minimum Spanning Tree Edges:");
        for (int i = 1; i < graph.length; i++) {
            System.out.println(parent[i] + " - " + i);
        }
    }

    private static int[] primMST(int[][] graph) {
        int n = graph.length;

        int[] parent = new int[n]; // Store the parent of each vertex in the MST
        int[] key = new int[n]; // Key values used to pick minimum weight edges
        boolean[] visited = new boolean[n]; // Track visited vertices

        Arrays.fill(key, INF); // Initialize key values with infinity
        Arrays.fill(visited, false); // Mark all vertices as not visited

        PriorityQueue<Node> pq = new PriorityQueue<>(Comparator.comparingInt(node -> node.key));
        pq.offer(new Node(0, 0)); // Start with the first vertex as the source
        key[0] = 0; // Set key value of source vertex to 0

        while (!pq.isEmpty()) {
            Node node = pq.poll();
            int u = node.vertex;

            visited[u] = true; // Mark the vertex as visited

            for (int v = 0; v < n; v++) {
                int weight = graph[u][v];

                if (!visited[v] && weight != INF && weight < key[v]) {
                    // Update key value and enqueue the neighboring vertex
                    pq.offer(new Node(v, weight));
                    parent[v] = u;
                    key[v] = weight;
                }
            }
        }

        return parent;
    }

    static class Node {
        int vertex;
        int key;

        Node(int vertex, int key) {
            this.vertex = vertex;
            this.key = key;
        }
    }
}





###### 7) Implement Greedy search algorithm for Kruskal's Minimal 
Spanning Tree Algorithm ######

#KruskalMST.java





import java.util.*;

public class KruskalMST {
    public static void main(String[] args) {
        // Create the graph
        int[][] graph = {
                {0, 2, 0, 6, 0},
                {2, 0, 3, 8, 5},
                {0, 3, 0, 0, 7},
                {6, 8, 0, 0, 9},
                {0, 5, 7, 9, 0}
        };

        // Find the minimum spanning tree
        List<Edge> mst = kruskalMST(graph);

        // Print the minimum spanning tree edges
        System.out.println("Minimum Spanning Tree Edges:");
        for (Edge edge : mst) {
            System.out.println(edge.src + " - " + edge.dest + ", Weight: " + edge.weight);
        }
    }

    private static List<Edge> kruskalMST(int[][] graph) {
        int n = graph.length;

        // Create a list to store the edges of the minimum spanning tree
        List<Edge> mst = new ArrayList<>();

        // Create a priority queue to store all the edges sorted by weight
        PriorityQueue<Edge> pq = new PriorityQueue<>(Comparator.comparingInt(edge -> edge.weight));

        // Add all edges to the priority queue
        for (int src = 0; src < n; src++) {
            for (int dest = src + 1; dest < n; dest++) {
                int weight = graph[src][dest];
                if (weight != 0) {
                    pq.offer(new Edge(src, dest, weight));
                }
            }
        }

        // Create a disjoint set to keep track of the connected components
        DisjointSet ds = new DisjointSet(n);

        while (!pq.isEmpty() && mst.size() < n - 1) {
            // Get the edge with the minimum weight from the priority queue
            Edge edge = pq.poll();
            int src = edge.src;
            int dest = edge.dest;

            // Check if adding this edge creates a cycle
            if (ds.find(src) != ds.find(dest)) {
                // Include the edge in the minimum spanning tree
                mst.add(edge);

                // Merge the two sets into one
                ds.union(src, dest);
            }
        }

        return mst;
    }

    static class Edge {
        int src;
        int dest;
        int weight;

        Edge(int src, int dest, int weight) {
            this.src = src;
            this.dest = dest;
            this.weight = weight;
        }
    }

    static class DisjointSet {
        int[] parent;
        int[] rank;

        DisjointSet(int n) {
            parent = new int[n];
            rank = new int[n];

            // Initialize parent and rank arrays
            for (int i = 0; i < n; i++) {
                parent[i] = i;
                rank[i] = 0;
            }
        }

        int find(int x) {
            if (parent[x] != x) {
                // Path compression: Make the parent of x the root of its subtree
                parent[x] = find(parent[x]);
            }
            return parent[x];
        }

        void union(int x, int y) {
            int xRoot = find(x);
            int yRoot = find(y);

            if (xRoot == yRoot) {
                return;
            }

            // Union by rank: Attach the smaller rank tree under the root of the higher rank tree
            if (rank[xRoot] < rank[yRoot]) {
                parent[xRoot] = yRoot;
            } else if (rank[xRoot] > rank[yRoot]) {
                parent[yRoot] = xRoot;
            } else {
                parent[yRoot] = xRoot;
                rank[xRoot]++;
            }
        }
    }
}







###### 8) Implement Greedy search algorithm for Dijkstra's Minimal
Spanning Tree Algorithm  ######


#DijkstraMST.java


import java.util.*;

public class DijkstraMST {
    private static final int INF = Integer.MAX_VALUE;

    public static void main(String[] args) {
        // Create the graph
        int[][] graph = {
                {0, 2, INF, 6, INF},
                {2, 0, 3, 8, 5},
                {INF, 3, 0, INF, 7},
                {6, 8, INF, 0, 9},
                {INF, 5, 7, 9, 0}
        };

        // Find the minimum spanning tree
        List<Edge> mst = dijkstraMST(graph);

        // Print the minimum spanning tree edges
        System.out.println("Minimum Spanning Tree Edges:");
        for (Edge edge : mst) {
            System.out.println(edge.src + " - " + edge.dest + ", Weight: " + edge.weight);
        }
    }

    private static List<Edge> dijkstraMST(int[][] graph) {
        int n = graph.length;

        List<Edge> mst = new ArrayList<>(); // Store the edges of the minimum spanning tree
        int[] parent = new int[n]; // Store the parent of each vertex in the MST
        int[] dist = new int[n]; // Store the cumulative weight to reach each vertex

        Arrays.fill(dist, INF); // Initialize all distances with infinity
        Arrays.fill(parent, -1); // Mark all vertices as unvisited

        dist[0] = 0; // Set the distance of the source vertex to 0

        PriorityQueue<Node> pq = new PriorityQueue<>(Comparator.comparingInt(node -> node.dist));
        pq.offer(new Node(0, 0)); // Add the source vertex to the priority queue

        while (!pq.isEmpty()) {
            Node node = pq.poll();
            int u = node.vertex;

            if (parent[u] != -1) {
                // Add the edge to the minimum spanning tree
                mst.add(new Edge(parent[u], u, graph[parent[u]][u]));
            }

            for (int v = 0; v < n; v++) {
                int weight = graph[u][v];

                if (weight != INF && weight < dist[v]) {
                    // Update the distance and parent of the neighboring vertex
                    dist[v] = weight;
                    parent[v] = u;
                    pq.offer(new Node(v, dist[v]));
                }
            }
        }

        return mst;
    }

    static class Node {
        int vertex;
        int dist;

        Node(int vertex, int dist) {
            this.vertex = vertex;
            this.dist = dist;
        }
    }

    static class Edge {
        int src;
        int dest;
        int weight;

        Edge(int src, int dest, int weight) {
            this.src = src;
            this.dest = dest;
            this.weight = weight;
        }
    }
}



###### 9) Implement Greedy search algorithm for Single-Source Shortest Path Problem  ######



#Dijkstra.java



import java.util.*;

public class Dijkstra {

    private static final int INF = Integer.MAX_VALUE;

    public static void main(String[] args) {
        // Example graph represented as an adjacency matrix
        int[][] graph = {
            {0, 4, 2, 0, 0, 0},
            {4, 0, 1, 5, 0, 0},
            {2, 1, 0, 8, 10, 0},
            {0, 5, 8, 0, 2, 6},
            {0, 0, 10, 2, 0, 3},
            {0, 0, 0, 6, 3, 0}
        };

        int source = 0; // Source vertex

        int[] shortestDistances = dijkstra(graph, source);

        // Print shortest distances from the source vertex to all other vertices
        System.out.println("Shortest Distances from source vertex " + source + ":");
        for (int i = 0; i < shortestDistances.length; i++) {
            System.out.println("Vertex " + i + ": " + shortestDistances[i]);
        }
    }

    private static int[] dijkstra(int[][] graph, int source) {
        int n = graph.length;
        int[] distances = new int[n]; // Shortest distances from the source vertex
        boolean[] visited = new boolean[n]; // Mark vertices as visited

        Arrays.fill(distances, INF); // Initialize distances with infinity
        distances[source] = 0; // Set distance of source vertex to 0

        for (int i = 0; i < n - 1; i++) {
            int minVertex = findMinVertex(distances, visited); // Find the vertex with the minimum distance

            visited[minVertex] = true; // Mark the vertex as visited

            for (int j = 0; j < n; j++) {
                // Update distances of adjacent vertices if a shorter path is found
                if (graph[minVertex][j] != 0 && !visited[j]
                        && distances[minVertex] + graph[minVertex][j] < distances[j]) {
                    distances[j] = distances[minVertex] + graph[minVertex][j];
                }
            }
        }

        return distances;
    }

    private static int findMinVertex(int[] distances, boolean[] visited) {
        int minVertex = -1;
        for (int i = 0; i < distances.length; i++) {
            if (!visited[i] && (minVertex == -1 || distances[i] < distances[minVertex])) {
                minVertex = i;
            }
        }
        return minVertex;
    }
}







