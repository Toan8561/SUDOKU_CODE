#pragma once

#include <iostream>
#include <stack>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "Definitions.h"

struct Node {
    Node* top;
    Node* bottom;
    Node* left;
    Node* right;
    Node* colHeader;
    int row;
    int column;
    int value;
    bool header;
    Node() {
        top = bottom = left = right = colHeader = NULL;
        row = column = value = -1;
        header = false;
    }
    Node(int r, int c, int v) {
    top = bottom = left = right = colHeader = NULL;
    row = r;
    column = c;
    value = v;
    header = false;
    }
};

class SudokuMatrix {
public:
    SudokuMatrix();
    ~SudokuMatrix();
    bool initialize();

    bool AddColumn(Node* newNode);

    void print();

    std::stack<Node>* solve(const char* filename);

private:
    Node* Root;
    std::stack<Node>* workingSolution; 
    bool Solved; 
    int totalCompetition;

    bool AddColumnHelp(Node* newNode, Node* r);

    bool isEmpty(); 

    void deleteMatrix();

    void cover(Node* r); 

    void uncover(Node* r); 

    Node* find(Node* find);

    bool solve();

    Node* chooseNextColumn(int& count);
};