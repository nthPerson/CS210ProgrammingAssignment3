#include <iostream>
#include <vector>

using namespace std;

class Edge {
public:
    int startVertex;
    int endVertex;
    int weight;

    // TODO Default constructor for Edge?
//    Edge() {
//        startVertex = -1;
//        endVertex = -1;
//        weight = 0;
//    }

    Edge(int weight, int startVertex, int endVertex) : weight(weight), startVertex(startVertex), endVertex(endVertex) {};
//    Edge(int weight, int startVertex, int endVertex) {
//        this->weight = weight;
//        this->startVertex = startVertex;
//        this->endVertex = endVertex;
//    }

    bool operator>(const Edge &otherEdge) const {
        // logic is reversed for use in min heap, smaller weight = higher priority
        return weight < otherEdge.weight;
//        return weight > otherEdge.weight;
    }

    bool operator<(const Edge &otherEdge) const {
        // logic is reversed for use in min heap, greater weight = lower priority
        return weight > otherEdge.weight;
//        return weight < otherEdge.weight;
    }

    bool operator==(const Edge &otherEdge) const {
        return weight == otherEdge.weight;
    }

    void print() {
        // TODO ? add endline in Edge::print()
        cout << startVertex << " -> " << endVertex << " : " << weight;
    }
};

// not sure if this will be necessary for Prim's
class Vertex {
public:
    int predecessor;
    int minWeight;
};

class Graph {
private:
    int numVertices;
    // 2D vector for adjacency matrix
    vector<vector<int>> adjacencyMatrix;

public:
    // make a new graph with a given number of vertices, all weights are initialized to 0
    Graph(int vertices) : numVertices(vertices), adjacencyMatrix(vertices, vector<int>(vertices, 0)) {}

    // make a new graph from any 2D vector
    Graph(const vector<vector<int>> &sampleGraph) {
        numVertices = sampleGraph.size();
        adjacencyMatrix = sampleGraph;
    }

//    Graph(int vertices, int sampleGraph[][5]) : numVertices(vertices), adjacencyMatrix(vertices, vector<int>(vertices, 0)) {
//        // convert sample graph to Graph type to allow for more flexibility for Prim's algo input
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
    T *top;
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
//    PriorityQueue() {
//        top = nullptr;
//
//
//    }

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

void primsMST(Graph &graph) {
    int numVertices = graph.getNumVertices();
    vector<bool> visited(numVertices, false);
    PriorityQueue<Edge> minEdgePQueue;
    vector<Edge> minSpanTree;
    int weight;

    // pick an arbitrary start vertex and add it to the lis of visited vertices
    int currentIndex = 0;
    visited[currentIndex] = true;

    // add neighbors of start vertex to priority queue
    for (int i = 0; i < numVertices; ++i) {
        if (graph.hasEdge(currentIndex, i)) {
            weight = graph.getEdgeWeight(currentIndex, i);
            minEdgePQueue.push(Edge(weight, currentIndex, i));
        }
    }

    // build the rest of the MST
    while (!minEdgePQueue.isEmpty() && minSpanTree.size() < numVertices - 1) {
        // minimum top will always be the top of the p-queue
        Edge minEdge = minEdgePQueue.getTop();
        minEdgePQueue.deleteTop();

        // TODO finish writing method

    }

}



int main() {
    vector<vector<int>> graph;

    int G[5][5] = { {0, 3, 65, 0, 0},
                    {3, 0, 85, 20, 45},
                    {65, 85, 0, 41, 77},
                    {0, 20, 41, 0, 51},
                    {0, 45, 77, 51, 0} };

    // convert sample adjacency matrix to vector to input into prims method

    // convert given 2D array graph to 2D vector to test new generalized Graph constructor
    vector<vector<int>> convertedSampleGraph(5, vector<int>(5));

    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            convertedSampleGraph[i][j] = G[i][j];
        }
    }

    // test PriorityQueue and Graph
    PriorityQueue<Edge> *pqTest;

//    Graph *testGraph = new Graph(5, G);
    Graph *testGraph = new Graph(convertedSampleGraph);
    Edge edge1(20, 1, 2);
    Edge edge2(20, 3, 4);
    Edge edge3(5, 5, 6);
    Edge edge4(10, 7, 8);

    testGraph->print();

//    if (pqTest->isEmpty()) {
//        cout << "PQ is empty." << endl;
//    } else {
//        cout << "isEmpty does not work." << endl;
//    }

    //test Edge::print()
//    cout << "Edge1: " << endl;
//    edge2.print();

    // test addEdge()
//    cout << "Graph before adding top:" << endl;
//    testGraph->print();
//
//    cout << "Graph after adding top:" << endl;
//    testGraph->addEdge(4, 1, 99999);
//    testGraph->print();

    // test get numVertices()
//    cout << "Number of vertices in testGraph: " << testGraph->getNumVertices() << endl;

    // test hasEdge
//    if (testGraph->hasEdge(0, 1)) {
//        cout << "Edge at (0, 1) confirmed." << endl;
//        cout << "Weight of top at (0, 1): " << testGraph->getEdgeWeight(0, 1) << endl;
//    }


    // Compare edge1 with edge2
//    if (edge1 > edge2) {
//        cout << "Edge1 is greater than Edge2" << endl;
//    } else if (edge1 == edge2) {
//        cout << "Edge1 is equal to Edge 2" << endl;
//    } else {
//        cout << "Edge1 is not greater than Edge2" << endl;
//    }
//
//    // Compare edge3 with edge1
//    if (edge3 < edge1) {
//        cout << "Edge3 is less than Edge1" << endl;
//    } else {
//        cout << "Edge3 is not less than Edge1" << endl;
//    }
//
//    // Comparing edges with equal weights
//    if (edge1 > edge4) {
//        cout << "Edge1 is greater than Edge4" << endl;
//    } else if (edge1 < edge4) {
//        cout << "Edge1 is less than Edge4" << endl;
//    } else {
//        cout << "Edge1 and Edge4 are equal in weight" << endl;
//    }


//    testGraph->print();



    return 0;
}
