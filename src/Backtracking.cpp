#include <iostream>
#include <fstream>
#include <string.h>
#include <queue>
#include <stdlib.h>
#include <time.h>
#include <list>
#include <chrono>
#include <ctime>
#include <stack>
#include <bits/stdc++.h>
#include <string>
#include <cstring>
#include "SudokuMatrix.h"
#include "Definitions.h"

#define SIZE 9
using namespace std;

#define MAX_DIGITS 10

int matrix[9][9] = {
    {0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0}
};


//function to print sudoku
void print_sudoku() {
    int i,j;
    for (i = 0; i < SIZE; i++) {
        cout << "|";
        for(j = 0; j < SIZE; j++) {
            cout <<" "<<matrix[i][j]<<" "<<"|";
        }
        cout << "\n+-----------+-----------+-----------+\n";
    }
}


//------------------Back tracking------------------//
//function to check if all cells are assigned or not
//if there is any unassigned cell then this function will change the values of row and col accordingly
bool number_unassigned(int *row, int *col) {
    for(int i = 0; i < SIZE; i++) {
        for(int j = 0; j < SIZE; j++) {
            if(matrix[i][j] == 0) {
                *row = i;
                *col = j;
                return true;
            }
        }
    }
    return false;
}

//function to check if we can put a value in a paticular cell or not
bool is_safe(int n, int r, int c) {
    for(int i = 0; i < SIZE; i++) {
        if(matrix[r][i] == n)
            return false;
    }
    for(int i = 0; i < SIZE; i++) {
        if(matrix[i][c] == n)
            return false;
    }
    int row_start = (r/3)*3;
    int col_start = (c/3)*3;
    for(int i = row_start; i < row_start + 3; i++) {
        for(int j = col_start; j < col_start + 3; j++) {
            if(matrix[i][j] == n)
                return false;
        }
    }
    return true;
}

//function to solve sudoku using backtracking
bool backtracking() {
    int row, col;
    if(!number_unassigned(&row, &col))
        return true;
    for(int i = 1; i <= SIZE; i++) {
        if(is_safe(i, row, col)) {
            matrix[row][col] = i;
            if(backtracking())
                return true;
            matrix[row][col]=0;
        }
    }
    return false;
}
//------------------------------------------------//

//------------------Bitwise------------------//
bool solve(int r, int c, int board[9][9], int submatrixDigits[3][3], int rowDigits[9], int columnDigits[9]) {
	if (r == 9) {
		return true;
	}
	if (c == 9) {
		return solve(r + 1, 0, board, submatrixDigits, rowDigits, columnDigits);
	}
	
	if (board[r][c] == 0) {
		for (int i = 1; i <= 9; i++) {
			int digit = 1 << (i - 1);
			if (!((submatrixDigits[r / 3][c / 3] & digit) || (rowDigits[r] & digit) || (columnDigits[c] & digit))) {
				submatrixDigits[r / 3][c / 3] |= digit;
				rowDigits[r] |= digit;
				columnDigits[c] |= digit;
				board[r][c] = i;
				if (solve(r, c + 1, board, submatrixDigits, rowDigits, columnDigits)) {
					return true;
				} else {
					submatrixDigits[r / 3][c / 3] &= ~digit;
					rowDigits[r] &= ~digit;
					columnDigits[c] &= ~digit;
					board[r][c] = 0;
				}
			}
		}
		return false;
	}
	return solve(r, c + 1, board, submatrixDigits, rowDigits, columnDigits);
}

bool SolveSudoku(int board[9][9]) {
	int submatrixDigits[3][3];
	int columnDigits[9];
	int rowDigits[9];

	for (int i = 0; i < 3; i++)
		memset(submatrixDigits[i], 0, 3 * sizeof(int));
	memset(rowDigits, 0, 9 * sizeof(int));
	memset(columnDigits, 0, 9 * sizeof(int));
	
	for (int i = 0; i < 9; i++)
		for (int j = 0; j < 9; j++)
			if (board[i][j] > 0) {
				int value = 1 << (board[i][j] - '1');
				submatrixDigits[i / 3][j / 3] |= value;
				rowDigits[i] |= value;
				columnDigits[j] |= value;
			}

	if (solve(0, 0, board, submatrixDigits, rowDigits, columnDigits))
		return true;
	else
		return false;
}

void printMatrix(int grid[9][9]) {
	for (int row = 0; row < 9; row++) {
	    for (int col = 0; col < 9; col++)
			cout << grid[row][col] << "\t";
		cout << "\n\n";
	}
}
//-------------------------------------------//

//----------------Exact cover----------------//
#include "SudokuMatrix.h"

SudokuMatrix::SudokuMatrix() {
    Root = new Node();
    Root->header=true;
    Root->right=Root->left=Root->top=Root->bottom=Root;
    workingSolution = new std::stack<Node>();
    Solved = false;
}

SudokuMatrix::~SudokuMatrix() {
    deleteMatrix(); 
    delete Root;
    delete workingSolution;
}

bool SudokuMatrix::AddColumn(Node* newNode) {
    if (!newNode->header)
        return false;
    else
        return AddColumnHelp(newNode,Root);
}

bool SudokuMatrix::AddColumnHelp(Node* newNode, Node* r) {
    if (r->right == Root && r != newNode) {
        r->right->left = newNode;
        newNode->right = r->right;
        newNode->left = r;
        r->right = newNode;
        return true;
    }
    else if (r == newNode)
        return false;
    else
        return AddColumnHelp(newNode,r->right);
}


void SudokuMatrix::print() {
    if (Root==NULL)
        return;

    int count = 0;
    int colCount = 0;
    Node* printNode = Root->right;
    Node* printNodeNext = Root->right;
    while(printNodeNext != Root) {
        colCount++;
        printNode = printNodeNext->bottom;
        while(!printNode->header) {
            printNode=printNode->bottom;
            count++;
        }
        printNodeNext = printNodeNext->right;
    }

    std::cout << "Counted " << count << " values, " << colCount << " Column Headers" << std::endl;
}

void SudokuMatrix::deleteMatrix() {
    Node* deleteNextCol = Root->right;
    Node* deleteNextRow;
    Node* temp;
    while(deleteNextCol != Root) {
        deleteNextRow = deleteNextCol->bottom;
        while(!deleteNextRow->header) {
            temp = deleteNextRow->bottom;
            delete deleteNextRow;
            deleteNextRow = temp;
        }
        temp = deleteNextCol->right;
        delete deleteNextCol;
        deleteNextCol = temp;
    }
}

std::stack<Node>* SudokuMatrix::solve(const char* filename) {
    Solved = false;
    totalCompetition = 0;
    while(!workingSolution->empty())
        workingSolution->pop();
    std::ifstream fin;
    fin.open(filename);
    if (fin.fail()) {
        std::cout << "Error, could not open " << filename << " for reading" << std::endl;
        fin.close();
        return NULL;
    }
    int nextVal;

    Node* toFind = NULL;
    Node* insertNext = NULL;
    std::stack<Node*> puzzleNodes;
    Node* colNode,*nextRowInCol,*rowNode;
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            fin >> nextVal;
            if (fin.fail() || nextVal > MATRIX_SIZE || nextVal < 0) {
	            std::cout << "Invalid Sudoku Puzzle specified" << std::endl;
	            fin.close();
	            return NULL;
            }

            if (nextVal != 0) {
	            toFind = new Node(i,j,nextVal-1);
	            insertNext = find(toFind);
	            if (insertNext == NULL) {
	                std::cout<<"Error in Sudoku Puzzle " << i << ", " << j << " val= " << nextVal << std::endl;
	                fin.close();
	                return NULL;
	            }
	            colNode = insertNext->colHeader;
	            nextRowInCol = insertNext;
	            cover(colNode);

	            rowNode = nextRowInCol->right;
	            while(rowNode != nextRowInCol) {
	                cover(rowNode->colHeader);
	                rowNode = rowNode->right;
                }

	            puzzleNodes.push(insertNext);
	            workingSolution->push(*insertNext);
	            delete toFind;
            }
        }
    }
    fin.close();

    std::cout << "Solving..." << std::endl;

    if (solve())
        std::cout << "Puzzle solved successfully!" << std::endl;
    else
        std::cout << "Puzzle not solveable!" << std::endl;

    while(!puzzleNodes.empty()) {
        colNode = (puzzleNodes.top())->colHeader;
        nextRowInCol = (puzzleNodes.top());

        rowNode = nextRowInCol->right;
        while(rowNode != nextRowInCol) {
            uncover(rowNode->colHeader);
            rowNode = rowNode->right;
        }
        uncover(colNode);
        puzzleNodes.pop();
    }

    std::stack<Node> temp,*toRet;
    while(!workingSolution->empty()) {
        temp.push(workingSolution->top());
        workingSolution->pop();
    }
    toRet = new std::stack<Node>();
    while(!temp.empty()) {
        toRet->push(temp.top());
        temp.pop();
    }

    return toRet;
}

void SudokuMatrix::cover(Node* r) {
    Node *RowNode, *RightNode, *ColNode = r->colHeader;
    ColNode->right->left = ColNode->left;
    ColNode->left->right = ColNode->right;
    for(RowNode = ColNode->bottom; RowNode != ColNode; RowNode = RowNode->bottom) {
        for(RightNode = RowNode->right; RightNode != RowNode; RightNode = RightNode->right) {
            RightNode->top->bottom = RightNode->bottom;
            RightNode->bottom->top = RightNode->top;
        }
    }
}

void SudokuMatrix::uncover(Node* r) {
    Node *RowNode, *LeftNode,*ColNode = r->colHeader;
    for(RowNode = ColNode->top; RowNode != ColNode; RowNode = RowNode->top) {
    for(LeftNode = RowNode->left; LeftNode!=RowNode; LeftNode = LeftNode->left) {
        LeftNode->top->bottom = LeftNode;
        LeftNode->bottom->top = LeftNode;
        }
    }
    ColNode->right->left = ColNode;
    ColNode->left->right = ColNode;
}

bool SudokuMatrix::isEmpty() {
    return (Root->right == Root);
}

bool SudokuMatrix::solve() {
    if (isEmpty())
        return true;

    int numCols;
    Node* nextCol = NULL;
    nextCol = chooseNextColumn(numCols);

    if (numCols < 1)
        return false;

    totalCompetition += numCols;

    Node* nextRowInCol = nextCol->bottom;
    Node* rowNode;
    cover(nextCol);

    while(nextRowInCol != nextCol && !Solved) {
        workingSolution->push(*nextRowInCol);
        rowNode = nextRowInCol->right;
        while(rowNode != nextRowInCol) {
            cover(rowNode->colHeader);
            rowNode = rowNode->right;
        }

        Solved=solve();
        if (!Solved) {
            workingSolution->pop();
        }

        rowNode = nextRowInCol->right;
        while(rowNode != nextRowInCol) {
            uncover(rowNode->colHeader);
            rowNode = rowNode->right;
        }

        nextRowInCol = nextRowInCol->bottom;
    }
    uncover(nextCol);
    return Solved;
}

Node* SudokuMatrix::chooseNextColumn(int& count) {
    Node* currentBest = Root->right;
    int best = -1;
    int tempCount = 0;

    Node* next = currentBest->bottom;
    Node* nextCol = currentBest;
    while(nextCol != Root) {
        next = nextCol->bottom;
        tempCount = 0;
        while(next != nextCol) { 
            if (next == next->bottom) {
	            std::cout<<"Err!" << std::endl;
            }
            tempCount++;
            next = next->bottom;
        }
        if (tempCount < best || best == -1) {
            currentBest = nextCol;
            best = tempCount;
        }
        nextCol = nextCol->right;
    }

    if (currentBest == Root) {
        std::cout << "Attempted to choose column from empty matrix!" << std::endl;
        exit(-1);
    }

    count = best;
    return currentBest;
}

Node* SudokuMatrix::find(Node* find) {
    Node* rightNode,*bottomNode;
    rightNode = Root->right;
    while(rightNode != Root) {
        bottomNode = rightNode->bottom;
        while(bottomNode != rightNode) {
            if (bottomNode->row == find->row && bottomNode->column == find->column && bottomNode->value == find->value) {
	            return bottomNode;
            }
            bottomNode = bottomNode->bottom;
        }
        rightNode = rightNode->right;
    }

    return NULL;
}

bool SudokuMatrix::initialize() {
    Node* matrix[MAX_ROWS][MAX_COLS];
    for (int i = 0; i < MAX_ROWS; i++) {
        for (int j = 0; j < MAX_COLS; j++) {
            matrix[i][j] = NULL;
        }
    }

    int row=0;
    Node *rowNode,*colNode,*cellNode,*boxNode;
    for (int i=0;i<MATRIX_SIZE;i++) {
        for (int j=0;j<MATRIX_SIZE;j++) {
            for (int k=0;k<MATRIX_SIZE;k++) {
	            row = (i*COL_OFFSET+j*MATRIX_SIZE+k);
	            rowNode=matrix[row][ROW_OFFSET+(i*MATRIX_SIZE+k)] = new Node(i,j,k);
	            colNode=matrix[row][COL_OFFSET+(j*MATRIX_SIZE+k)] = new Node(i,j,k);
	            cellNode=matrix[row][CELL_OFFSET+(i*MATRIX_SIZE+j)] = new Node(i,j,k);
	            boxNode=matrix[row][BOX_OFFSET+((i/ROW_BOX_DIVISOR + j/COL_BOX_DIVISOR * COL_BOX_DIVISOR)*MATRIX_SIZE+k)] = new
	            Node(i,j,k);
	            rowNode->right =colNode;
	            rowNode->left=boxNode;
	            colNode->left=rowNode;
	            colNode->right=cellNode;
	            cellNode->left=colNode;
	            cellNode->right=boxNode;
	            boxNode->left=cellNode;
	            boxNode->right=rowNode;
            }
        }
    }
    Node* nextColHeader;
    Node* nextColRow;
    for (int j = 0; j < MAX_COLS; j++) {
        nextColHeader = new Node();
        nextColHeader->header=true;
        nextColHeader->top=nextColHeader->bottom=nextColHeader->left=nextColHeader->right=nextColHeader->colHeader=nextColHeader;
        nextColRow = nextColHeader;

        for (int i=0;i<MAX_ROWS;i++) {
            if (matrix[i][j] != NULL) {
	            matrix[i][j]->top=nextColRow;
	            nextColRow->bottom = matrix[i][j];
	            matrix[i][j]->bottom = nextColHeader;
	            nextColHeader->top = matrix[i][j];
	            matrix[i][j]->colHeader = nextColHeader;
	            nextColRow = matrix[i][j];
            }
        }
        if (nextColHeader->bottom == nextColHeader) {
            std::cout<<"Err! column has no rows! col:" << j << std::endl;
        }
        if(!AddColumn(nextColHeader)) {
            std::cout << "Error in adding column to matrix" << std::endl;
            for (int i=0;i<MAX_ROWS;i++) {
	            for (int j=0;j<MAX_COLS;j++) {
	                if (matrix[i][j] != NULL)
	                    delete matrix[i][j];
	            }
            }
            return false;
        }
    }
    return true;
}

bool initialize(SudokuMatrix* &m);

void getinfo(void){
    fstream myFile;
        myFile.open("log1.txt", ios::in);
        if (myFile.is_open()){
            string line;
            while (getline(myFile, line)){
                char p[line.length()];
                p[0] = line[0];
                if(p[0] == '*'){
                    for (int i = 1; i < sizeof(p); i++) {
                     p[i] = line[i];
                     cout << p[i];
                    }
                    cout << endl;
                }
                else break;
            }
        }  
        myFile.close();   
}

int getprob(string line){
    int n;
    cout << "enter id number" << endl;
    cin >> n;
    char name[]="+PROBLEM_NO.";
    char num_char[MAX_DIGITS + sizeof(char)];
    sprintf(num_char, "%d", n);
    string names = string(name);
    string num_chars = string(num_char);
    string temp= names+num_chars;
    int m = temp.length();
    char id[m+1];
    strcpy(id, temp.c_str());
    int a = line.length();
    char char_array[a + 1];
    strcpy(char_array, line.c_str());
    int result = strcmp(id, char_array);
    return result;
}

void getmatrix(void){
    int i=0;
    fstream myFile;
    myFile.open("log1.txt", ios::in);
    if (myFile.is_open()){
            string line;
            while (getline(myFile, line)){
                char p[line.length()];
                p[0] = line[0];
                if(p[0] == '*' || p[0] == '%' || p[0] == '$' || p[0] == '=') continue;
                else if(p[0] == '+') {
                    if(getprob(line)==0){
                        cout << "Found the log" << endl;
                        continue;
                    }
                    else {
                        cout << "Not found int log" << endl;
                        continue;
                    }
                }
                else if(p[0] == '#'){
                    for (int i = 1; i < sizeof(p); i++) {
                     p[i] = line[i];
                     cout << p[i];
                    }
                    cout << " is loading" << endl;
                    continue;
                    }
                else if(p[0] == '>'){
                    cout << "getting data" <<endl;
                    int m=0;
                    for (int j = 1; j < sizeof(p);) {
                     p[j] = line[j];
                     if(p[j]!=' ') {
                        matrix[i][m]=p[j]-'0';
                        m++;
                        j++;
                        }
                     else j++;
                    }
                    i++;
                    }
                else break;
            }
        }  
        myFile.close();   
}

void write_solution(int type){
    int i,j;
    fstream myFile;
    myFile.open("log1.txt", ios::app);
    if (myFile.is_open()){
        myFile << "\n";
        for (i = 0; i < SIZE; i++) {
            myFile << "<";
            for(j = 0; j < SIZE; j++) {
                myFile <<matrix[i][j]<<" ";
            }
            myFile << "\n";
        }
        switch (type) {
            case 1: {
                myFile << "Solution: Backtracking \n" ;
                break;
                }
            case 2: {
                myFile << "Solution: Bitwise \n" ;
                break;
            }
            default: break;
        }
        myFile <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n\n";
        cout << "would you like to run again"<<endl;
        cout << "1 -> Yes    0 -> No"<<endl;
        cout << "Enter: ";
        int o;
        cin >> o;
        if(o==1) myFile <<"=SOLUTION \n";
    }
    myFile.close(); 
}

//-------------------------------------------//

int main(void) {
    getinfo();
    getmatrix();
    srand(time(0));
    cout << "Solving the sudoku:\n";
    SudokuMatrix* m = new SudokuMatrix();
    string userInput;
    bool userContinues = true;
    stack<Node>* solution;
    Node next;

    if (!(m->initialize())) {
        cout << "Could not initialize matrix" << endl;
        delete m;
        return 1;
    }
    cout << "Initialized matrix" << endl << endl;
     cout << "+-----------+-----------+-----------+\n";
    int puzzle[MATRIX_SIZE][MATRIX_SIZE]; 
    string puzzleToSolve = "example.txt";
    // cout << "\nFilename of puzzle to solve: ";
    print_sudoku();
    cout << "We have these algorithm:\n";
    cout << "1. Backtracking\n";
    cout << "2. Bitwise\n";
    cout << "3. Exact cover\n";;
    cout << "Choose your algorithm (1-3): ";
    int n;
    cin >> n;
    switch (n) {
    case 1:
        if (backtracking()) {
            cout << "Solution: " << endl;
            cout << "\n+-----------+-----------+-----------+\n";
            print_sudoku();
            auto end = chrono::system_clock::now();
            write_solution(n);
        }
        else
            cout << "No solution exists\n";
        break;
    
    case 2:
        if (SolveSudoku(matrix)){
            cout << "Solution: " << endl;
            cout << "\n+-----------+-----------+-----------+\n";
            print_sudoku();
            auto end = chrono::system_clock::now();
            write_solution(n);
        }
        else
            cout << "No solution exists\n";
        break;

    case 3:{
        fstream myFile;
        myFile.open("log1.txt", ios::app);
        if (myFile.is_open()){
        myFile << "\n";
        solution = m->solve(puzzleToSolve.c_str());
        if (solution != NULL) {
            while(!solution->empty()) {
	            next = solution->top();
	            puzzle[next.row][next.column] = next.value+1;
	            solution->pop();
            }
        } 
        else cout << "Solution could not be found" << endl;
        cout << "Solution: " << endl;
        cout << "\n+-----------+-----------+-----------+\n";
        
        for (int i = 0; i < MATRIX_SIZE; i++) {
            cout <<"|";
            myFile << ">";
            for (int j = 0; j < MATRIX_SIZE; j++) {
	            cout <<" " << puzzle[i][j] << " " << "|";
                myFile << puzzle[i][j]<<" ";
            }
             cout << "\n+-----------+-----------+-----------+\n";
            myFile << "\n";
        }
        myFile << "Solution: Exact cover\n" ;
        myFile <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n\n";
        cout << "would you like to run again"<<endl;
        cout << "1 -> Yes    0 -> No"<<endl;
        cout << "Enter: ";
        int o;
        cin >> o;
        if(o==1) myFile <<"=SOLUTION \n";

        delete solution;
        solution = NULL;
        break;
        }
        myFile.close();
    }

    default:
        break;
    }

    if (solution != NULL) {
        while(!solution->empty())
            solution->pop();
        delete solution;
    }
    if (m != NULL) 
        delete m;
    return 0;
}