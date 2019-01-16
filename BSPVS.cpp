// BSPVS.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

//#include "pch.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>

//Declare I/O variables
std::ofstream os1;
std::ofstream os2;
std::ifstream is;
const std::string openFileName("tp");
//Declare Global variables
int totalNodeNum = 0;
int totalEdgeNum = 0;
int limitTravelTime = 0;
int startTime = 0;
int resultHappinessCredit1 = 0;
int resultHappinessCredit2 = 0;
int nodeOfTree = 0;
//Declare Constant
const int empty = -1;
const int nodeLimit = 1290000;
//Declare Class
class Node {
public:
	Node(std::string n, int h, int o, int c, int i, int* t) : name(n), happy(h), open(o), close(c), index(i), table(t) {};
	Node(std::string n, int h, int o, int c, int i) : name(n), happy(h), open(o), close(c), index(i), table(NULL) {};
	Node(std::string n, int h, int o, int c) : name(n), happy(h), open(o), close(c), index(-1), table(NULL) {};
	Node() : name(""), happy(-1), open(-1), close(-1), index(-1), table(NULL) {};

	const std::string name; const int happy; const int open; const int close; const int index; int* table; std::vector<Node*> neighbor;
};
class Edge {
public:
	Edge(std::string n1, std::string n2, int w) :name1(n1), name2(n2), weight(w) {};
	Edge() :name1(""), name2(""), weight(-1) {};
	const std::string name1; const std::string name2; const int weight;
};
class GraphM : public Node, public Edge {
public:
	GraphM() : nodeNum(0), edgeNum(0) {};
	
	Node node(std::string name);
	Node node(int i);
	Node* nodeP(std::string name);
	Node* nodeP(int i);
	Edge edge(std::string name1, std::string name2);
	Edge edge(unsigned int i1, unsigned int i2);
	int distance(std::string name1, std::string name2);
	int distance(unsigned int i1, unsigned int i2);
	bool addNode(Node nd);
	bool addEdge(Edge ed);
	void show();
	std::vector<Node> nodeList; int nodeNum; int edgeNum;
};
//Global Class Variables
GraphM graph;

class NodeA {
public:
	NodeA() : index(empty), node(NULL), accumHapp(empty), accumDis(empty), parent(NULL), child(0), everPass(NULL) {};
	NodeA(Node* nd, int aH, int aD, NodeA* p) : index(nd->index), node(nd), accumHapp(aH), accumDis(aD), parent(p), child(0), everPass(NULL) {};
	NodeA(Node* nd, int aH, int aD, NodeA* p, bool* eP) : index(nd->index), node(nd), accumHapp(aH), accumDis(aD), parent(p), child(0), everPass(eP) {};
	NodeA(Node* nd) : index(nd->index), node(nd), accumHapp(nd->happy), accumDis(0), parent(NULL), child(0), everPass() {
		everPass = new bool[totalNodeNum] {false};
		everPass[nd->index] = true;
		for (int i = 0; i < totalNodeNum; i++) {
			if (everPass[i] == true) { printf("1"); }
			else if (everPass[i] == false) { printf("0"); }
		}
	};
	
	NodeA* newChild(Node* nd, int disLimit);
	NodeA* qDeleteExcept(std::queue<NodeA*>* q);
	void del();
	void contiDel(std::queue<NodeA*>* leaf);
	int accumHaITime();
	int index; Node* node; int accumHapp; int accumDis; NodeA* parent; int child; bool* everPass;
};
class CmpElem {
public:
	CmpElem(int idx) : index(idx), cmpIdx(0) {
		Node* nd = graph.nodeP(idx);
		unsigned int n = nd->neighbor.size();
		int weightSum = 0;
		for (unsigned int i = 0; i < n; i++) {
			weightSum = weightSum + graph.distance(nd->neighbor.at(i)->index, idx);
		}
		cmpIdx = (float)(nd->happy)*(float)(n) / (float)(weightSum);
	}
	int index; float cmpIdx;
};
//Constant
const Node nullNode;
const Edge nullEdge;
const NodeA nullNodeA;
//Declare Utility Functions
bool openFileTemp();
void readFile();
Node getNode(std::string n, int h, int o, int c) { Node nd(n, h, o, c); return nd; }
Edge getEdge(std::string n1, std::string n2, int h) { Edge ed(n1, n2, h); return ed; }
float divisor(int i) { if (i == 0) { return 0.5; } return i; }
void isCorrect1(NodeA* end, int finalHI, int finalDis);
void isCorrect2(NodeA* end, int finalHI, int finalDis);

//Algorithm1
//BFS + Trim(Maybe)
class trimRatio {
public:

};

std::pair<NodeA*, NodeA*> BFS(int idx, int dis) {
	printf("Node %d Start: \n", idx);
	NodeA* start = new NodeA(graph.nodeP(idx));
	std::queue<NodeA*> q;
	//std::queue<NodeA*> qLeaf;
	std::queue<NodeA*> qRecord;
	std::vector<NodeA*> vecEnd;
	//std::vector<NodeA*> vecEndTime;
	q.push(start);
	//printf("Push Node: %d\n", start->index);

	while (!q.empty()) {
		NodeA* now = q.front();
		std::vector<NodeA*> vecTrim;
		std::vector<Node*> vecNeighbor(now->node->neighbor);
		float avgPerform = 0;
		int notNull = 0;
		bool hasChild = false;

		for (unsigned int i = 0; i < vecNeighbor.size(); i++) {
			if (nodeOfTree > nodeLimit) { break; }
			Node* newNode{ vecNeighbor.at(i) };
			NodeA* newNodeA = now->newChild(newNode, limitTravelTime);
			if (newNodeA != NULL) {
				//printf("New Child Node: %d OF %d AHI: %d, AD: %d\n", now->node->neighbor.at(i)->index, now->index, newNodeA->accumHapp, newNodeA->accumDis);
				vecTrim.push_back(newNodeA);
				avgPerform = avgPerform + (float)(newNodeA->accumHapp) / divisor(newNodeA->accumDis);
				notNull++;
				hasChild = true;
			}
		}
		//if (!hasChild) { qLeaf.push(now); }
		avgPerform = avgPerform / divisor(notNull);
		//Trim
		float trimStandard = avgPerform;
		int trimmedNum = 0;
		//printf("Avg: %f, Trim Standard: %f\n", avgPerform, trimStandard);
		for (unsigned int i = 0; i < vecTrim.size(); i++) {
			NodeA* temp = vecTrim.at(i);
			float performP = (float)(temp->accumHapp) / divisor(temp->accumDis);
			if (performP < trimStandard) {
				//printf("Trimmed: %f\n", performP);
				if (temp->node->open < startTime + temp->accumDis && temp->node->close > startTime + temp->accumDis) {
					nodeOfTree++;
					if (temp->index == idx) { vecEnd.push_back(temp);/*vecEndTime.push_back(temp);*/ }
				}
				else {

					delete(temp);
					vecTrim.erase(vecTrim.begin() + i); trimmedNum++;//Optimize
				}
			}
			else {
				q.push(temp);
				//qRecord.push(temp);
				//printf("Push Node: %d\n", temp->index);
				nodeOfTree++;
				//std::cout << nodeOfTree << "\n";
				if (temp->index == idx) { vecEnd.push_back(temp); }
			}
		}
		float trimRatio = (float)(trimmedNum) / divisor(vecTrim.size());
		//printf("Trim Ratio: %d/%d = %f", trimmedNum, vecTrim.size(), trimRatio);
		q.pop();
	}
	//Pick Best Node
	NodeA* MaxP;//No Time Limit
	if (vecEnd.size() > 0) {
		MaxP = vecEnd.front();
		for (unsigned int i = 1; i < vecEnd.size(); i++) {
			if (vecEnd.at(i)->accumHapp > MaxP->accumHapp) { MaxP = vecEnd.at(i); }
		}
	}
	else {
		MaxP = start;
	}
	NodeA* MaxPTime;//Time Limit
	int haITime = 0;
	if (vecEnd.size() > 0) {
		MaxPTime = vecEnd.front();
		haITime = vecEnd.front()->accumHaITime();
		for (unsigned int i = 1; i < vecEnd.size(); i++) {
			NodeA* temp = vecEnd.at(i);
			int tempHaITime = vecEnd.at(i)->accumHaITime();
			if (tempHaITime > haITime) { MaxPTime = temp; haITime = tempHaITime; }
		}
	}
	else {
		haITime = start->accumHaITime();
		MaxPTime = start;
	}
	NodeA* MaxPTimeAlt = new NodeA;
	MaxPTimeAlt->index = MaxPTime->index;
	MaxPTimeAlt->node = MaxPTime->node;
	MaxPTimeAlt->accumHapp = haITime;
	MaxPTimeAlt->accumDis = MaxPTime->accumDis;
	MaxPTimeAlt->parent = MaxPTime->parent;
	MaxPTimeAlt->child = MaxPTime->child;
	MaxPTimeAlt->everPass = MaxPTime->everPass;

	float r = (float)(MaxP->node->happy);
	int td = 0;
	int neighborOfMaxP = MaxP->node->neighbor.size();
	for (int i = 0; i < neighborOfMaxP; i++) {
		td = td + graph.distance(MaxP->node->neighbor.at(i)->index, MaxP->index);
	}
	r = r / (float)(td) * (float)(neighborOfMaxP);
	std::cout << "\n" << nodeOfTree << "\n";
	printf("MaxP: index=%d, R(HI:%d ,EN:%d, AVGEW:%f)=%f, accumHapp=%d, accumDis=%d\n", MaxP->index, MaxP->node->happy, neighborOfMaxP, (float)(td) / (float)(neighborOfMaxP), r, MaxP->accumHapp, MaxP->accumDis);
	//MaxP->contiDel(&qLeaf);
	//MaxP->qDeleteExcept(&qRecord);
	//nodeOfTree = 0;
	
	return std::make_pair(MaxP, MaxPTimeAlt);
}

void plan() {
	std::vector<CmpElem> vecCmp;
	for (int i = 0; i < totalNodeNum; i++) {
		CmpElem e(i);
		vecCmp.push_back(e);
	}
	std::sort(vecCmp.begin(), vecCmp.end(), [](CmpElem a, CmpElem b) {return a.cmpIdx > b.cmpIdx; });
	std::pair<NodeA*, NodeA*> maxNode = BFS(vecCmp.at(0).index, limitTravelTime);
	for (int i = 1; i < totalNodeNum; i++) {
		std::pair<NodeA*, NodeA*> temp = BFS(vecCmp.at(i).index, limitTravelTime);		NodeA* tempF = temp.first;
		NodeA* tempS = temp.second;
		if (maxNode.first->accumHapp < temp.first->accumHapp) { maxNode.first = temp.first; }
		if (maxNode.second->accumHapp < temp.second->accumHapp) { maxNode.second = temp.second; }
	}
	//Output ans1
	std::cout << "\n" << maxNode.first->accumHapp << " " << maxNode.first->accumDis << "\n";
	os1 << maxNode.first->accumHapp << " " << maxNode.first->accumDis << "\n";
	std::stack<NodeA*> stck;
	NodeA* temp = maxNode.first;
	while (1) {
		stck.push(temp);
		if (temp->parent == NULL) { break; }
		temp = temp->parent;
	}
	printf("%d output: \n", stck.size());

	unsigned int stckSize = stck.size();
	for (unsigned int i = 0; i < stckSize; i++) {
		NodeA* tempNodeA = stck.top();
		int accumTime = tempNodeA->accumDis + startTime;
		std::cout << i << " ";
		std::cout << tempNodeA->node->name << " " << accumTime << " " << accumTime << "\n";
		os1 << tempNodeA->node->name << " " << accumTime << " " << accumTime << "\n";

		stck.pop();
	}
	isCorrect1(maxNode.first, maxNode.first->accumHapp, maxNode.first->accumDis);

	//Output ans2
	std::cout << "\n" << maxNode.second->accumHapp << " " << maxNode.second->accumDis << "\n";
	os2 << maxNode.second->accumHapp << " " << maxNode.second->accumDis << "\n";
	std::stack<NodeA*> stck2;
	NodeA* temp2 = maxNode.second;
	while (1) {
		stck2.push(temp2);
		if (temp2->parent == NULL) { break; }
		temp2 = temp2->parent;
	}
	printf("%d output: \n", stck2.size());

	unsigned int stckSize2 = stck2.size();
	for (unsigned int i = 0; i < stckSize2; i++) {
		NodeA* tempNodeA = stck2.top();
		int accumTime = tempNodeA->accumDis + startTime;
		std::cout << i << " ";
		std::cout << tempNodeA->node->name << " " << accumTime << " " << accumTime << "\n";
		os2 << tempNodeA->node->name << " " << accumTime << " " << accumTime << "\n";

		stck2.pop();
	}
	isCorrect2(maxNode.second, maxNode.second->accumHapp, maxNode.second->accumDis);
}

int main(int argc, char** argv) {
	/*
	std::string p(argv[1]);
	  std::string path = "./" + p + "/tp.data";
	  is.open(path);
	  path = "";
	  path = "./" + p + "/ans1.txt";
	  os1.open(path);
	path = "";
	  path = "./" + p + "/ans2.txt";
	  os2.open(path);*/
	if (!openFileTemp()) { std::cout << "Fail to open file" << '\n'; return 0; }
	readFile();
	plan();

	is.close();
	os1.close();
	os2.close();
}

//Implement
//Utility Implement--------------------------------------------------------------------------------------
bool openFileTemp() {
	is.open(openFileName + ".data");
	os1.open("ans1.txt");
	os2.open("ans2.txt");
	return is.is_open() && os1.is_open() && os2.is_open();
}
void readFile() {
	//Graph Abstract
	is >> totalNodeNum; is >> totalEdgeNum; is >> limitTravelTime; is >> startTime;
	std::cout << totalNodeNum << " " << totalEdgeNum << " " << limitTravelTime << " " << startTime << '\n';

	//Graph Nodes Detail
	for (int i = 0; i < totalNodeNum; i++) {
		std::string name; int happyIndex; int bussiTimeStar; int bussiTimeEnd;
		is >> name; is >> happyIndex; is >> bussiTimeStar; is >> bussiTimeEnd;
		std::cout << name << " " << happyIndex << " " << bussiTimeStar << " " << bussiTimeEnd << '\n';
		graph.addNode(getNode(name, happyIndex, bussiTimeStar, bussiTimeEnd));
	}
	//Graph Edges Detail
	for (int i = 0; i < totalEdgeNum; i++) {
		std::string name1; std::string name2; int weight;
		is >> name1; is >> name2; is >> weight;
		std::cout << name1 << " " << name2 << " " << weight << '\n';
		graph.addEdge(getEdge(name1, name2, weight));
	}
}
void isCorrect1(NodeA * end, int finalHI, int finalDis){
	NodeA* temp;
	int hi = 0;
	int dis = 0;
	bool* ep = new bool[totalNodeNum];
	if (end != NULL) { temp = end; }
	else { return; }
	for (int i = 0; i < totalNodeNum; i++) {
		ep[i] = false;
	}
	while (temp != NULL) {
		if (ep[temp->index] == false) { hi = hi + temp->node->happy; ep[temp->index] = true; }
		if(temp->parent != NULL){
			dis = dis + graph.distance(temp->parent->index, temp->index);
			temp = temp->parent;
		}
		else {
			break;
		}
	}
	if (finalDis != dis) { printf("Incorrect Dis, correct: %d\n", dis); }
	if (finalHI != hi) { printf("Incorrect HI, correct: %d\n", hi); }
	if (dis > limitTravelTime) { printf("Out Of Bound, dis: %d\n", dis); }
}
void isCorrect2(NodeA * end, int finalHI, int finalDis){
	NodeA* temp;
	int hi = 0;
	int dis = 0;
	bool* ep = new bool[totalNodeNum];
	if (end != NULL) { temp = end; }
	else { return; }
	for (int i = 0; i < totalNodeNum; i++) {
		ep[i] = false;
	}
	while (temp != NULL) {
		int current = startTime + temp->accumDis;
		if (ep[temp->index] == false && temp->node->open < current && current < temp->node->close) { hi = hi + temp->node->happy; ep[temp->index] = true; }
		if (temp->parent != NULL) {
			dis = dis + graph.distance(temp->parent->index, temp->index);
			temp = temp->parent;
		}
		else {
			break;
		}
	}
	if (finalDis != dis) { printf("Incorrect Dis, correct: %d\n", dis); }
	if (finalHI != hi) { printf("Incorrect HI, correct: %d\n", hi); }
	if (dis > limitTravelTime) { printf("Out Of Bound, dis: %d\n", dis); }
}
//Utility Implement^--------------------------------------------------------------------------------------
//GraphM Implement---------------------------------------------------------------------------------
Node GraphM::node(std::string name) {
	for (unsigned int i = 0; i < this->nodeList.size(); i++) {
		if (this->nodeList.at(i).name == name) { return this->nodeList.at(i); }
	}
	return nullNode;
}
Node GraphM::node(int i) {
	unsigned int j = i;
	if (j < nodeList.size()) { return nodeList.at(i); }
	else { return nullNode; }
}
Node* GraphM::nodeP(std::string name) {
	for (unsigned int i = 0; i < this->nodeList.size(); i++) {
		if (nodeList.at(i).name == name) { return &(nodeList.at(i)); }
	}
	return new Node(nullNode);
}
Node* GraphM::nodeP(int i) {
	unsigned int j = i;
	if (j < nodeList.size()) { return &(nodeList.at(i)); }
	else { return new Node(nullNode); }
}
Edge GraphM::edge(std::string name1, std::string name2) {
	int idx1 = -1, idx2 = -1; bool n1Found = false, n2Found = false;
	if (name1 == name2) { return nullEdge; }
	for (unsigned int i = 0; i < nodeList.size(); i++) {
		if (nodeList.at(i).name == name1) {
			idx1 = i; n1Found = true;
			if (n1Found && n2Found) { return getEdge(name2, name1, nodeList.at(idx1).table[idx2]); }
		}
		else if (nodeList.at(i).name == name2) {
			idx2 = i; n2Found = true;
			if (n1Found && n2Found) { return getEdge(name1, name1, nodeList.at(idx2).table[idx1]); }
		}
	}
	return nullEdge;
}
Edge GraphM::edge(unsigned int i1, unsigned int i2) {
	if (i1 < nodeList.size() && i2 < nodeList.size() && i1 != i2) {
		if (i2 > i1) { return getEdge(nodeList.at(i1).name, nodeList.at(i2).name, nodeList.at(i2).table[i1]); }
		else if (i1 > i2) { return getEdge(nodeList.at(i2).name, nodeList.at(i1).name, nodeList.at(i1).table[i2]); }
		else { return getEdge(nodeList.at(i1).name, nodeList.at(i2).name, nodeList.at(i2).table[i1]); }
		//else{return nullEdge;}
	}
	else {
		return nullEdge;
	}
}
int GraphM::distance(std::string name1, std::string name2) {
	int idx1 = -1, idx2 = -1; bool n1Found = false, n2Found = false;
	if (name1 == name2) { return -1; }
	for (unsigned int i = 0; i < this->nodeList.size(); i++) {
		if (this->nodeList.at(i).name == name1) {
			idx1 = i; n1Found = true;
			if (n1Found && n2Found) { return this->nodeList.at(idx1).table[idx2]; }
		}
		else if (this->nodeList.at(i).name == name2) {
			idx2 = i; n2Found = true;
			if (n1Found && n2Found) { return this->nodeList.at(idx2).table[idx1]; }
		}
	}
	return -1;
}
int GraphM::distance(unsigned int i1, unsigned int i2) {
	if (i1 < this->nodeList.size() && i2 < this->nodeList.size() && i1 != i2) {
		if (i2 > i1) { return this->nodeList.at(i2).table[i1]; }
		else if (i1 > i2) { return this->nodeList.at(i1).table[i2]; }
		else { return this->nodeList.at(i2).table[i1]; }
		//else{return nullEdge;}
	}
	else {
		return -1;
	}
}
bool GraphM::addNode(Node nd) {
	nodeNum++;
	Node nds(nd.name, nd.happy, nd.open, nd.close, this->nodeList.size(), new int[this->nodeNum + 1]);
	for (int i = 0; i < nodeNum + 1; i++) { nds.table[i] = -1; } this->nodeList.push_back(nds);
	return true;
}
bool GraphM::addEdge(Edge ed) {
	int idx1 = -1, idx2 = -1; bool n1Found = false, n2Found = false;
	if (ed.name1 == ed.name2) { return false; }
	for (unsigned int i = 0; i < this->nodeList.size(); i++) {
		if (this->nodeList.at(i).name == ed.name1) {
			idx1 = i; n1Found = true;
			if (n1Found && n2Found) {
				this->nodeList.at(idx1).neighbor.push_back(&(this->nodeList.at(idx2)));
				this->nodeList.at(idx2).neighbor.push_back(&(this->nodeList.at(idx1)));
				edgeNum++; this->nodeList.at(idx1).table[idx2] = ed.weight;
				return true;
			}
		}
		else if (this->nodeList.at(i).name == ed.name2) {
			idx2 = i; n2Found = true;
			if (n1Found && n2Found) {
				this->nodeList.at(idx1).neighbor.push_back(&(this->nodeList.at(idx2)));
				this->nodeList.at(idx2).neighbor.push_back(&(this->nodeList.at(idx1)));
				edgeNum++; this->nodeList.at(idx2).table[idx1] = ed.weight;
				return true;
			}
		}
	}
	return false;
}
void GraphM::show() {
	for (int i = 0; i < nodeNum; i++) {
		for (int j = 0; j < nodeNum; j++) {
			std::cout << graph.edge(i, j).weight << " ";
		}
		std::cout << '\n';
	}
}
//GraphM Implement^---------------------------------------------------------------------------------
NodeA* NodeA::newChild(Node* nd, int disLimit) {
	bool* everPassN = new bool[graph.nodeNum];
	for (int i = 0; i < graph.nodeNum; i++) { everPassN[i] = this->everPass[i]; }
	int nowDis = this->accumDis + graph.distance(nd->index, this->index);

	if (nowDis <= disLimit) {
		this->child++;
		if (!this->everPass[nd->index]) {
			everPassN[nd->index] = true;
			return new NodeA(nd, this->accumHapp + nd->happy, nowDis, this, everPassN);
		}
		else {
			return new NodeA(nd, this->accumHapp, nowDis, this, everPassN);
		}
	}
	else { return NULL; }
}
void NodeA::del() {
	this->parent->child--;
	delete(this);
}
NodeA* NodeA::qDeleteExcept(std::queue<NodeA*>* q) {
	NodeA* temp(this);
	while (temp->parent != NULL) {
		this->child = -1;
		temp = temp->parent;
	}
	q->pop();
	q->pop();
	while (!q->empty()) {
		if (q->front()->child != -1) {
			delete (q->front());
		}
		q->pop();
	}
	return this;
}
void NodeA::contiDel(std::queue<NodeA*>* leaf) {//Enter a Queue of All Leaves

	//NodeA* parentOfFnt = leaf->front()->parent;
	//NodeA* fnt = leaf->front();
	NodeA* parentOfFnt;
	NodeA* fnt;
	
	while (!leaf->empty()) {
		//std::cout << "Not Empty\n";
		fnt = leaf->front();
		parentOfFnt = leaf->front()->parent;
		if (fnt != this) {
			if (fnt->child == 0) { 
				//std::cout << "fnt" << fnt->index << " Child = 0\n";
				parentOfFnt->child--; delete(fnt);
			}
			if (parentOfFnt->child == 0) {
				//std::cout << "parent child = 0\n";
				leaf->push(parentOfFnt); 
			}
			//std::cout << "Leaf" << leaf->front()->index << " Pop " << leaf->size() << "\n";
			leaf->pop();
		}
		else { 
			//std::cout << "Exception Leaf Pop\n";
			leaf->pop();
		}
		/*if (leaf->size() > 0) {
			fnt = leaf->front();
			parentOfFnt = leaf->front()->parent;
		}*/
	}
}

int NodeA::accumHaITime() {
	NodeA* temp = this;
	int HaI = 0;
	bool* everTravel;
	everTravel = new bool[totalNodeNum];
	for (int i = 0; i < totalNodeNum; i++) {
		everTravel[i] = false;
	}
	while (temp != NULL) {
		int timeNow = startTime + temp->accumDis;
		if (temp->node->open < timeNow && timeNow < temp->node->close && everTravel[temp->index] == false) {
			HaI = HaI + temp->node->happy;
			everTravel[temp->index] = true;
		}
		if (temp->parent == NULL) { break; }
		
		temp = temp->parent;
	}
	return HaI;
}


// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
