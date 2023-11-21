//
// Created by rober on 11/21/2023.
//

#include "testMethods.h"
// PriorityQueue test methods

void testPush() {
    PriorityQueue<int> pq;
    pq.push(10);
    pq.push(5);
    pq.push(15);
//    PriorityQueue<Edge> pq;
//    pq.push(Edge(10, 0, 1));
//    pq.push(Edge(5, 1, 2));
//    pq.push(Edge(15, 2, 3));

    assert(pq.getTop() == 5);  // The smallest element should be on top
    cout << "Test Push: Passed" << endl;
}
void testDeleteTop() {
    PriorityQueue<int> pq;
    pq.push(10);
    pq.push(5);
    pq.push(15);
    pq.deleteTop();  // Removes 5

    assert(pq.getTop() == 10);  // Next smallest element should be on top
    cout << "Test Delete Top: Passed" << endl;
}
void testHeapProperty() {
    PriorityQueue<int> pq;
    pq.push(20);
    pq.push(30);
    pq.push(10);
    pq.push(5);
    pq.push(15);

    assert(pq.getTop() == 5);
    pq.deleteTop();
    assert(pq.getTop() == 10);
    pq.deleteTop();
    assert(pq.getTop() == 15);

    cout << "Test Heap Property: Passed" << endl;
}
void testEmptyQueue() {
    PriorityQueue<int> pq;
    try {
        pq.deleteTop();  // Should handle empty queue
    } catch (const std::exception& e) {
        cout << "Test Empty Delete: Passed" << endl;
    }

    try {
        pq.getTop();  // Should handle empty queue
    } catch (const std::exception& e) {
        cout << "Test Empty Get Top: Passed" << endl;
    }
}
// Edge type tests
void testPushEdges() {
    PriorityQueue<Edge> pq;
    pq.push(Edge(10, 0, 1));
    pq.push(Edge(5, 1, 2));
    pq.push(Edge(15, 2, 3));

    assert(pq.getTop().weight == 5);  // Edge with smallest weight should be on top
    cout << "Test Push Edges: Passed" << endl;
}
void testDeleteTopEdges() {
    PriorityQueue<Edge> pq;
    pq.push(Edge(10, 0, 1));
    pq.push(Edge(5, 1, 2));
    pq.push(Edge(15, 2, 3));
    pq.deleteTop();  // Removes edge with weight 5

    assert(pq.getTop().weight == 10);  // Next smallest weight should be on top
    cout << "Test Delete Top Edges: Passed" << endl;
}void testHeapPropertyEdges() {
    PriorityQueue<Edge> pq;
    pq.push(Edge(20, 0, 1));
    pq.push(Edge(30, 1, 2));
    pq.push(Edge(10, 2, 3));
    pq.push(Edge(5, 3, 4));
    pq.push(Edge(15, 4, 5));

    assert(pq.getTop().weight == 5);
    pq.deleteTop();
    assert(pq.getTop().weight == 10);
    pq.deleteTop();
    assert(pq.getTop().weight == 15);

    cout << "Test Heap Property Edges: Passed" << endl;
}
void testEmptyQueueEdges() {
    PriorityQueue<Edge> pq;
    try {
        pq.deleteTop();  // Should handle empty queue
    } catch (const std::exception& e) {
        cout << "Test Empty Delete Edges: Passed" << endl;
    }

    try {
        pq.getTop();  // Should handle empty queue
    } catch (const std::exception& e) {
        cout << "Test Empty Get Top Edges: Passed" << endl;
    }
}


// from main()


//    // test PriorityQueue and Graph
//    PriorityQueue<Edge> *pqTest;
//
//    Graph *testGraph = new Graph(5, G);
//    Graph *testGraph = new Graph(convertedSampleGraph);
//    Edge edge1(20, 1, 2);
//    Edge edge2(20, 3, 4);
//    Edge edge3(5, 5, 6);
//    Edge edge4(10, 7, 8);
//    testGraph->print();

//    // test PriorityQueue using integer type
//    testPush();
//    testDeleteTop();
//    testHeapProperty();
//    testEmptyQueue();
//    // PQ tests using Edge type
//    testPushEdges();
//    testDeleteTopEdges();
//    testHeapPropertyEdges();
//    testEmptyQueueEdges();


//    if (pqTest->isEmpty()) {
//        cout << "PQ is empty." << endl;
//    } else {
//        cout << "isEmpty does not work." << endl;
//    }

// test getVertexEdges
//    cout << "All edges for vertex 1:" << endl;
//    vector<Edge> v1Edges = testGraph->getVertexEdges(1);
//    printEdges(v1Edges);




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