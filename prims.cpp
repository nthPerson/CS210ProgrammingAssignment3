#include <iostream>
#include <vector>
#include <cassert>

using namespace std;

class Edge {
public:
    int fromVertex;
    int toVertex;
    int weight;

    // initialization list (constructor)
    Edge(int weight, int startVertex, int endVertex) : weight(weight), fromVertex(startVertex), toVertex(endVertex) {};

    // overloaded operators to facilitate comparison of Edges in priority queue
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
        cout << fromVertex << " - " << toVertex << " -> " << weight;
    }


};

class Graph {
private:
    int numVertices;
    // 2D vector adjacency matrix
    vector<vector<int>> adjacencyMatrix;

public:
    // make a new graph with a given number of vertices, all weights are initialized to 0
    explicit Graph(int vertices) : numVertices(vertices), adjacencyMatrix(vertices, vector<int>(vertices, 0)) {}

    // make a new graph from any 2D vector
    explicit Graph(const vector<vector<int>> &twoDimVector) {
        adjacencyMatrix = twoDimVector;
    }

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

    // used to find edges of currentVertex in primsMST()
    vector<Edge> getVertexEdges(int vertex) {
        vector<Edge> edges;
        for (int i = 0; i < adjacencyMatrix.size(); ++i) {
            if(adjacencyMatrix[vertex][i] != 0 && adjacencyMatrix[vertex][i] != INT_MAX) {
                edges.emplace_back(adjacencyMatrix[vertex][i], vertex, i);
            }
        }
        return edges;
    }

    void print() {
        // adjacencyMatrix.size() returns number of rows, not number of elements
        for (int i = 0; i < adjacencyMatrix.size(); ++i) {
            for (int j = 0; j < adjacencyMatrix.size(); ++j) {
                cout << adjacencyMatrix[i][j] << ", ";
            }
            cout << endl;
        }
    }

};

// used to maintain discovered edges and find minimum weight edge
template<typename T>
class PriorityQueue {
private:
    vector<T> heap;

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
    // vector to keep track of which vertices have been visited
    vector<bool> visited(numVertices, false);
    // priority queue to maintain discovered edges and produce minimum weight edge
    PriorityQueue<Edge> minEdgePQ;
    // MST of the given graph
    vector<Edge> minimumSpanningTree;
    int weight;

    // pick an arbitrary start vertex and mark it visited
    int currentVertex = 0;
    visited[currentVertex] = true;

    // add edges from start vertex to priority queue
    // getVertexEdges uses the Edge constructor to assign the fromVertex and toVertex properties of each edge
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
     * Step 1: Uncomment these three lines:
     *         int profMatrix[ENTER ARRAY][SIZE HERE] = { PASTE MATRIX HERE }
     *         Graph professorManjuGraph = createGraphFromMatrix(profMatrix)
     *         primsMST(professorManjuGraph)
     * Step 2: Copy/paste your 2D array into the indicated location and enter the size of the array in the declaration
     * Step 3: Run the program and check the console for output
     * step 4: Have a great day!
     * */

    // EXAMPLE USAGE
    int G[5][5] = { {0, 3, 65, 0, 0},
                    {3, 0, 85, 20, 45},
                    {65, 85, 0, 41, 77},
                    {0, 20, 41, 0, 51},
                    {0, 45, 77, 51, 0} };
    Graph g = createGraphFromMatrix(G);
    primsMST(g);

    int exampleAdjMatrix[6][6] = {{0,  5,  0,  23, 0,  0},
                                  {5,  0,  10, 41, 12, 30},
                                  {0,  10, 0,  0,  8,  0},
                                  {23, 41, 0,  0,  16, 14},
                                  {0,  12, 8,  16, 0,  27},
                                  {0,  30, 0,  14, 27, 0}};
    Graph exampleGraph = createGraphFromMatrix(exampleAdjMatrix);
    primsMST(exampleGraph);

    int bigAdjMatrix[9][9] = { {0, 6, 0, 0, 14, 2, 0, 0, 0},
                               {6, 0, 17, 0, 0, 12, 0, 0, 0},
                               {0, 17, 0, 60, 0, 13, 41, 0, 0},
                               {0, 0, 60, 0, 0, 0, 0, 0, 0},
                               {14, 0, 0, 0, 0, 0, 0, 0, 4},
                               {2, 12,13, 0, 0, 0, 0, 21, 0},
                               {0, 0, 41, 0, 0, 0, 0, 17, 0},
                               {0, 0, 0, 0, 0, 21, 17, 0, 0, },
                               {0, 0, 0, 0, 4, 0, 0, 0, 0}};
    Graph bigGraph = createGraphFromMatrix(bigAdjMatrix);
    primsMST(bigGraph);

    // ********** User's 2D array can go here ***************************************************************
    // ********** Uncomment lines, enter your matrix and its size in the indicated locations, then run ******
//    int profMatrix[ENTER ARRAY][SIZE HERE] = { PASTE 2D ARRAY/MATRIX HERE };
//
//    Graph professorManjuGraph = createGraphFromMatrix(profMatrix);
//    primsMST(professorManjuGraph);


    return 0;
}
