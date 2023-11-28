#include <iostream>
#include <vector>
#include <cassert>

using namespace std;

class Edge {
public:
//    int fromVertex = INT_MAX;
//    int toVertex = INT_MAX;
    int fromVertex;
    int toVertex;
    int weight;

    // Initialization list (constructor)
    Edge(int weight, int startVertex, int endVertex) : weight(weight), fromVertex(startVertex), toVertex(endVertex) {};
//    Edge(int weight, int fromVertex, int toVertex) {
//        this->weight = weight;
//        this->fromVertex = fromVertex;
//        this->toVertex = toVertex;
//    }

    bool operator>(const Edge &otherEdge) const {
        return weight > otherEdge.weight;
    }

    bool operator<(const Edge &otherEdge) const {
        return weight < otherEdge.weight;
    }

    bool operator==(const Edge &otherEdge) const {
        return weight == otherEdge.weight;
    }

    void print() {
//        cout << fromVertex << " -> " << toVertex << " " << weight;
        cout << fromVertex << " - " << toVertex << " -> " << weight;
    }


};

// not sure if this will be necessary for Prim's
class Vertex {
public:
    int vertexNumber;
//    string name;
//    int minWeight;

    Vertex(int vertexNumber): vertexNumber(vertexNumber) {};

    int getVertexNumber() {
        return vertexNumber;
    }
    void print() {
        cout << vertexNumber << endl;
    }
};

class Graph {
private:
    int numVertices;
    // 2D vector for adjacency matrix
    vector<vector<int>> adjacencyMatrix;

public:
    // make a new graph with a given number of vertices, all weights are initialized to 0
    Graph(int vertices) : numVertices(vertices), adjacencyMatrix(vertices, vector<int>(vertices, 0)) {}

    // make a new graph fromVertex any 2D vector
    explicit Graph(const vector<vector<int>> &twoDimVector) {
        numVertices = twoDimVector.size();
        adjacencyMatrix = twoDimVector;
    }


//    Graph(int vertices, int sampleGraph[][5]) : numVertices(vertices), adjacencyMatrix(vertices, vector<int>(vertices, 0)) {
//        // convert sample graph toVertex Graph type toVertex allow for more flexibility for Prim's algo input
//        for (int i = 0; i < 5; ++i) {
//            for (int j = 0; j < 5; ++j) {
//                this->adjacencyMatrix[i][j] = sampleGraph[i][j];
//            }
//        }
//    }

    void addEdge(int startVertex, int endVertex, int weight) {
        adjacencyMatrix[startVertex][endVertex] = weight;
        // since the graph is undirected, add reverse of top
        adjacencyMatrix[endVertex][startVertex] = weight;
    }

    int getNumVertices() {
        return adjacencyMatrix.size();
    }

    bool hasEdge(int row, int col) {
        if (adjacencyMatrix[row][col] > 0) {
            return true;
        }
    }

    int getEdgeWeight(int row, int col) {
        return adjacencyMatrix[row][col];
    }

    vector<Edge> getVertexEdges(int vertex) {
        vector<Edge> edges;
        for (int i = 0; i < adjacencyMatrix.size(); ++i) {
            if(adjacencyMatrix[vertex][i] != 0 && adjacencyMatrix[vertex][i] != INT_MAX) {
                edges.push_back(Edge(adjacencyMatrix[vertex][i], vertex, i));
            }
        }
        return edges;
    }

    void print() {
        // Graph.size() returns number of rows, not number of elements
        for (int i = 0; i < adjacencyMatrix.size(); ++i) {
            for (int j = 0; j < adjacencyMatrix.size(); ++j) {
                cout << adjacencyMatrix[i][j] << ", ";
            }
            cout << endl;
        }
    }

};

template<typename T>
class PriorityQueue {
private:
//    T top = heap[0];
    vector<T> heap;


//    ~PriorityQueue() {
//
//    }

    void downHeap(int index) {
        int leftChildIndex = 2 * index + 1;
        int rightChildIndex = 2 * index + 2;
        int smallestIndex = index;

        if (leftChildIndex < heap.size() && heap[leftChildIndex] < heap[smallestIndex]) {
            smallestIndex = leftChildIndex;
        }

        if (rightChildIndex < heap.size() && heap[rightChildIndex] < heap[smallestIndex]) {
            smallestIndex = rightChildIndex;
        }

        if (smallestIndex != index) {
            // swap
            swap(heap[index], heap[smallestIndex]);
            downHeap(smallestIndex);
        }
    }

    void upHeap(int index) {
        int parentIndex = (index - 1) / 2;

        if (index && heap[parentIndex] > heap[index]) {
            // swap
            swap(heap[index], heap[parentIndex]);
            upHeap(parentIndex);
        }
    }

public:

    bool isEmpty() {
        return heap.empty();
    }

    T getTop() {
        if (isEmpty()) {
//            cout << "Priority queue is empty, nothing to return" << endl;
            throw out_of_range("Priority queue is empty,");
        }
        return heap.front();
    }

    void push(const T &value) {
        heap.push_back(value);
        // add to end of array
        int index = heap.size() - 1;
        upHeap(index);
    }

    void deleteTop() {
        if (isEmpty()) {
//            cout << "Priority queue is empty" <<  endl;
            throw out_of_range("Priority queue is empty,");
        }
        heap[0] = heap.back();
        heap.pop_back();
        downHeap(0);
    }

    void print() {

    }

};

// original solution
void primsMST(Graph &graph) {
    int numVertices = graph.getNumVertices();
    // vector toVertex keep track of which vertices have been visited
    vector<bool> visited(numVertices, false);
    // priority queue toVertex maintain discovered edges and produce minimum weight edge
    PriorityQueue<Edge> minEdgePQ;
    // MST of the given graph
    vector<Edge> minimumSpanningTree;
    int weight;

    // pick an arbitrary start vertex and mark it visited
    int currentVertex = 0;
    visited[currentVertex] = true;

    // add edges from start vertex to priority queue
    // getVertexEdges uses the Edge constructor to assign the fromVertex and toVertex characteristics of each edge
    vector<Edge> edges = graph.getVertexEdges(currentVertex);
    for (Edge edge: edges) {
        minEdgePQ.push(edge);
    }


    // build the rest of the MST
    // while the priority queue isn't empty and there are less than V-1 vertices in the MST...
    while (minimumSpanningTree.size() < numVertices - 1) {
        // PQ is a min-heap so top element will always be the minimum weight edge

        // get next minimum weight edge and remove it from priority queue
        Edge minEdge = minEdgePQ.getTop();
        minEdgePQ.deleteTop();

        // check if the toVertex of current min edge has been visited,
        // skip adding this edge to MST if visited with 'continue' to avoid creating cycles
        int toVertex = minEdge.toVertex;
        if (visited[toVertex]) {
            continue;
        }

        // add edge to the MST and mark the toVertex as visited
        minimumSpanningTree.push_back(minEdge);
        // visited[toVertex] needs to be true
        visited[toVertex] = true;

        // add neighbors of newly visited vertex (toVertex) to priority queue
        vector<Edge> newEdges = graph.getVertexEdges(toVertex);
        for (Edge edge: newEdges) {
            minEdgePQ.push(edge);
        }

    }

    cout << "Minimum spanning tree:" << endl;
    cout << "(Start Vertex - End Vertex -> Weight)" << endl;
    for (Edge edge: minimumSpanningTree) {
        edge.print();
        cout << endl;
    }
}


void printEdges(const vector<Edge> &edges) {
    for (const Edge &edge: edges) {
        cout << edge.fromVertex << " - " << edge.toVertex << " -> " << edge.weight << endl;
    }
}

template<size_t ROWS, size_t COLS>
Graph createGraphFromMatrix(int (&twoDimensionalArray)[ROWS][COLS]) {
    // construct empty graph using the integer type constructor
    Graph graph(ROWS);

    // add edges to graph using the information stored in the adjacency matrix
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            // only add elements of adj matrix that are nonzero, cuz those are the edges
            if (i != j && twoDimensionalArray[i][j] != 0) {
                graph.addEdge(i, j, twoDimensionalArray[i][j]);
            }
        }
    }

    return graph;
}

int main() {

    /* INSTRUCTIONS FOR USE:
     * Step 1: Copy/paste the 2D array into indicated location (or wherever you want, you know how this works)
     * Step 2: Pass your array as an argument to createGraphFromMatrix() method
     * Step 3: Uncomment these two lines:
     *         Graph professorManjuGraph = createGraphFromMatrix()
     *         primsMST(professorManjuGraph)
     * Step 4: Run the program and check the console for output
     * step 5: Have a great day!
     * */

    // Your 2D array can go here, if you'd like


//    Graph professorManjuGraph = createGraphFromMatrix(/*YOUR 2D ARRAY GOES HERE*/);
//    primsMST(professorManjuGraph);


    int G[5][5] = { {0, 3, 65, 0, 0},
                    {3, 0, 85, 20, 45},
                    {65, 85, 0, 41, 77},
                    {0, 20, 41, 0, 51},
                    {0, 45, 77, 51, 0} };
    Graph graph1 = createGraphFromMatrix(G);
    primsMST(graph1);



//    int G[5][5] = { {0, 3, 65, 0, 0},
//                    {3, 0, 85, 20, 45},
//                    {65, 85, 0, 41, 77},
//                    {0, 20, 41, 0, 51},
//                    {0, 45, 77, 51, 0} };


    // convert given 2D array graph to 2D vector to construct Graph
//    vector<vector<int>> *sampleMatrix(5, vector<int>(5));
    vector<vector<int>> sampleMatrix(5, vector<int>(5));


    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            sampleMatrix[i][j] = G[i][j];
        }
    }
    auto *graph = new Graph(sampleMatrix);

    cout << "Sample matrix output:" << endl;
    primsMST(*graph);

    // different adjacency matrix
    int D[6][6] = { {0, 5, 0, 23, 0, 0},
                    {5, 0, 10, 41, 12, 30},
                    {0, 10, 0, 0, 8, 0},
                    {23, 41, 0, 0, 16, 14},
                    {0, 12, 8, 16, 0, 27},
                    {0, 30, 0, 14, 27, 0} };
    // convert to vector for construction of graph
    vector<vector<int>> differentMatrix(6, vector<int>(6));
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            differentMatrix[i][j] = D[i][j];
        }
    }
    auto *differentGraph = new Graph(differentMatrix);

    cout << "Different matrix output:" << endl;
    primsMST(*differentGraph);

    // another test matrix
    int bigAdjMatrix[9][9] = { {0, 6, 0, 0, 14, 2, 0, 0, 0},
                            {6, 0, 17, 0, 0, 12, 0, 0, 0},
                            {0, 17, 0, 60, 0, 13, 41, 0, 0},
                            {0, 0, 60, 0, 0, 0, 0, 0, 0},
                            {14, 0, 0, 0, 0, 0, 0, 0, 4},
                            {2, 12,13, 0, 0, 0, 0, 21, 0},
                            {0, 0, 41, 0, 0, 0, 0, 17, 0},
                            {0, 0, 0, 0, 0, 21, 17, 0, 0, },
                            {0, 0, 0, 0, 4, 0, 0, 0, 0}};
    vector<vector<int>> bigVectorMatrix(9, vector<int>(9));
    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < 9; ++j) {
            bigVectorMatrix[i][j] = bigAdjMatrix[i][j];
        }
    }
    auto *bigGraph = new Graph(bigVectorMatrix);

    cout << "Big graph output:" << endl;
    primsMST(*bigGraph);




    return 0;
}
