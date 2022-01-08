// application.cpp <Starter Code>
// Arthor: Hua Xu
//
// University of Illinois at Chicago
// CS 251: Spring 2021
// Project #7 - Openstreet Maps
// This Project is used for Navigating the uic campus and get the path
// and shortest distances
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//

#include <iostream>
#include <iomanip>  /*setprecision*/
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <queue>
#include <algorithm>
#include <stack>
#include "graph.h"
#include "tinyxml2.h"
#include "dist.h"
#include "osm.h"

using namespace std;
using namespace tinyxml2;


class prioritize {
public:
    bool operator() (const pair<long long, double>& p1,
        const pair<long long, double> p2) const {
        return p1.second > p2.second;
    }
};

// this function to read all the node in Footways and return a vector
vector<long long> readTheID(vector<FootwayInfo> Footways) {
    vector<long long> result;
    // loop through all the vector
    for (size_t i = 0; i < Footways.size(); i++) {
        for (size_t j = 0; j < Footways[i].Nodes.size(); j++) {
            result.push_back(Footways[i].Nodes[j]);
        }
    }
    return result;
}

// this function is Dijkstra shortest path algorithm implement by priority_queue
// and return a map with predecessor
map<long long, long long>  Dijkstra(graph< long long, double>& G,
    long long startV, vector<FootwayInfo> Footways,
    map<long long, double>& distances) {
    double INF = numeric_limits<double>::max();
    map<long long, long long> prev;
    priority_queue<pair<long long, double>, vector<pair<long long, double>>,
        prioritize> unvisitedQueue;
    set<long long> adjV;  // to store neighbors
    long long currentV;
    vector<long long> inGraph = readTheID(Footways);
    double edgeWeight, alternativePathDistance;
    for (size_t i = 0; i < inGraph.size(); i++) {
        distances[inGraph[i]] = INF;   // all distances is INF
        prev[inGraph[i]] = -10;     // all predecessor is a negtive number
        unvisitedQueue.push(make_pair(inGraph[i], INF));
    }
    unvisitedQueue.push(make_pair(startV, 0));
    distances[startV] = 0;
    while (!unvisitedQueue.empty()) {  // loop through all the priority_queue
        currentV = unvisitedQueue.top().first;
        unvisitedQueue.pop();
        adjV = G.neighbors(currentV);  // get the neighbors
        for (auto x : adjV) {   //loop through whole neighbors
            if (G.getWeight(currentV, x, edgeWeight)) {
                alternativePathDistance = distances[currentV] + edgeWeight;
                if (alternativePathDistance < distances[x]) {
                    distances[x] = alternativePathDistance;
                    unvisitedQueue.push(make_pair(x, alternativePathDistance));
                    prev[x] = currentV;
                }
            }
        }
    }
    return prev;
}

// this functiion search the building ID with abbreviation or part name
// return a pos of building in vector
int searchTheBuilding(string buildname, vector<BuildingInfo> Buildings) {
    for (size_t i = 0; i < Buildings.size(); i++) {
        if (Buildings[i].Abbrev == buildname) {
            return i;
        }
    }
    for (size_t i = 0; i < Buildings.size(); i++) {
        if (Buildings[i].Fullname.find(buildname) != string::npos) {
            return i;
        }
    }
    return -1;
}

// this fucntion help us find the Nearest building and return a ID.
// it uses lenaer search to find the cloest building
long long findNearest(BuildingInfo building, vector<FootwayInfo> Footways,
    map<long long, Coordinates> Nodes) {
    long long cloest;
    double lat, _long;
    lat = building.Coords.Lat;
    _long = building.Coords.Lon;
    cloest = Footways[0].Nodes[0];
    // calculate the distance between Footways and buildings
    double minDistance = distBetween2Points(lat, _long,
        Nodes[Footways[0].Nodes[0]].Lat,
        Nodes[Footways[0].Nodes[0]].Lon);
    double distance;
    // find the min;
    for (size_t i = 0; i < Footways.size(); i++) {
        for (size_t j = 0; j < Footways[i].Nodes.size(); j++) {
            distance = distBetween2Points(lat, _long,
                Nodes[Footways[i].Nodes[j]].Lat,
                Nodes[Footways[i].Nodes[j]].Lon);
            if (minDistance > distance) {
                minDistance = distance;
                cloest = Footways[i].Nodes[j];
            }
        }
    }
    return cloest;
}

// this functiion to add vertex based on Nodes.
void addVerticeG(map<long long, Coordinates> Nodes,
    graph< long long, double>& G) {
    for (auto i : Nodes) {
        if (!G.addVertex(i.first))
            cout << "**Error: unable to add vertex '" << i.first << "', why not?" << endl;
    }
}

// this function add the edge and calculate two node's distance
// we added A to B and also B to A
void addEdge(graph< long long, double>& G, vector<FootwayInfo> Footways,
    map<long long, Coordinates>  Nodes) {
    double lat1, long1, lat2, long2;
    double weight;
    for (size_t i = 0; i < Footways.size(); i++) {
        for (size_t j = 0; j < Footways[i].Nodes.size() - 1; j++) {
            lat1 = Nodes[Footways[i].Nodes[j]].Lat;  // get the latitude
            long1 = Nodes[Footways[i].Nodes[j]].Lon;
            lat2 = Nodes[Footways[i].Nodes[j + 1]].Lat;
            long2 = Nodes[Footways[i].Nodes[j + 1]].Lon;
            // calculate the distBetween2Points
            weight = distBetween2Points(lat1, long1, lat2, long2);
            if (!G.addEdge(Footways[i].Nodes[j], Footways[i].Nodes[j + 1], weight)) {
                cout << "**Error: unable to add edge (" << Footways[i].Nodes[j] << ",";
                cout << Footways[i].Nodes[j + 1] << "," << weight << "), why not?" << endl;
            }
            if (!G.addEdge(Footways[i].Nodes[j + 1], Footways[i].Nodes[j], weight)) {
                cout << "**Error: unable to add edge (" << Footways[i].Nodes[j] << ",";
                cout << Footways[i].Nodes[j + 1] << "," << weight << "), why not?" << endl;
            }
        }
    }
}

// this fucntion to do simple output work for navagation guidence
void pinrtGuidence(int result, int result2, long long& startingID,
    long long& endingID, vector<BuildingInfo>Buildings,
    vector<FootwayInfo> Footways, map<long long, Coordinates> Nodes) {
    startingID = findNearest(Buildings[result], Footways, Nodes);
    endingID = findNearest(Buildings[result2], Footways, Nodes);
    cout << "Starting point: " << endl;
    cout << " " << Buildings[result].Fullname << endl;
    cout << " (" << Buildings[result].Coords.Lat << ", ";
    cout << Buildings[result].Coords.Lon << ")" << endl;
    cout << "Destination point:" << endl;
    cout << " " << Buildings[result2].Fullname << endl;
    cout << " (" << Buildings[result2].Coords.Lat << ", ";
    cout << Buildings[result2].Coords.Lon << ")" << endl;
    cout << endl;
    cout << "Nearest start node:" << endl;
    cout << " " << startingID << endl;
    cout << " (" << Nodes[startingID].Lat << ", ";
    cout << Nodes[startingID].Lon << ")" << endl;
    cout << "Nearest destination node:" << endl;
    cout << " " << endingID << endl;
    cout << " (" << Nodes[endingID].Lat << ", ";
    cout << Nodes[endingID].Lon << ")" << endl;
    cout << endl;
    cout << "Navigating with Dijkstra..." << endl;
}

// generate shortest path map and distances, then print the
// shortest distances in miles and the path.
void printshortestPath(long long startingID, long long endingID,
    vector<FootwayInfo> Footways, graph< long long, double>& G) {
    map<long long, double> distances;
    map<long long, long long> prev;
    vector<long long> di;
    double INF = numeric_limits<double>::max();
    prev = Dijkstra(G, startingID, Footways, distances);
    if (distances[endingID] == INF) {
        cout << "Sorry, destination unreachable" << endl;
    }
    else {
        cout << "Distance to dest: " << distances[endingID] << " miles" << endl;
        long long IDs;
        stack<long long> mystack;
        IDs = endingID;
        cout << "Path: ";
        // put the path into the stack
        while (IDs != startingID) {
            mystack.push(IDs);
            IDs = prev[IDs];
        }
        cout << startingID;
        // output the stack
        while (!mystack.empty()) {
            cout << "->" << mystack.top();
            mystack.pop();
        }
        cout << endl;
    }
}

// the navagation engine, input the name of start and destination and output
// the path or other information
void shortestPath(vector<BuildingInfo>Buildings,
    vector<FootwayInfo> Footways, map<long long, Coordinates> Nodes,
    graph< long long, double>& G) {
    string startBuilding, destBuilding;
    cout << "Enter start (partial name or abbreviation), or #> ";
    getline(cin, startBuilding);
    while (startBuilding != "#") {
        cout << "Enter destination (partial name or abbreviation)> ";
        getline(cin, destBuilding);
        int result, result2;
        result = searchTheBuilding(startBuilding, Buildings);
        result2 = searchTheBuilding(destBuilding, Buildings);
        // we find both addresses are valid
        if (result != -1 && result2 != -1) {
            long long startingID, endingID;
            pinrtGuidence(result, result2, startingID, endingID, Buildings,
                Footways, Nodes);
            printshortestPath(startingID, endingID, Footways, G);
            // if start address is not valid
        }
        else if (result == -1 && result2 != -1) {
            cout << "Start building not found" << endl;
            // if destination is not valid
        }
        else if (result != -1 && result2 == -1) {
            cout << "Destination building not found" << endl;
        }
        else {
            cout << "Start building not found" << endl;
        }
        cout << endl;
        cout << "Enter start (partial name or abbreviation), or #> ";
        getline(cin, startBuilding);
    }
}


int main() {
    // maps a Node ID to it's coordinates (lat, lon)
    map<long long, Coordinates>  Nodes;
    // info about each footway, in no particular order
    vector<FootwayInfo>          Footways;
    // info about each building, in no particular order
    vector<BuildingInfo>         Buildings;
    XMLDocument                  xmldoc;

    cout << "** Navigating UIC open street map **" << endl;
    cout << endl;
    cout << std::setprecision(8);

    string def_filename = "map.osm";
    string filename;

    cout << "Enter map filename> ";
    getline(cin, filename);

    if (filename == "") {
        filename = def_filename;
    }

    //
    // Load XML-based map file
    //
    if (!LoadOpenStreetMap(filename, xmldoc)) {
        cout << "**Error: unable to load open street map." << endl;
        cout << endl;
        return 0;
    }

    //
    // Read the nodes, which are the various known positions on the map:
    //
    int nodeCount = ReadMapNodes(xmldoc, Nodes);

    //
    // Read the footways, which are the walking paths:
    //
    int footwayCount = ReadFootways(xmldoc, Footways);

    //
    // Read the university buildings:
    //
    int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

    //
    // Stats
    //

    assert(nodeCount == (int)Nodes.size());
    assert(footwayCount == (int)Footways.size());
    assert(buildingCount == (int)Buildings.size());
    cout << endl;
    cout << "# of nodes: " << Nodes.size() << endl;
    cout << "# of footways: " << Footways.size() << endl;
    cout << "# of buildings: " << Buildings.size() << endl;
    graph< long long, double> G;
    addVerticeG(Nodes, G);
    addEdge(G, Footways, Nodes);
    cout << "# of vertices: " << G.NumVertices() << endl;
    cout << "# of edges: " << G.NumEdges() << endl;
    cout << endl;
    shortestPath(Buildings, Footways, Nodes, G);
    //
    // done:
    //
    cout << "** Done **" << endl;
    return 0;
}
