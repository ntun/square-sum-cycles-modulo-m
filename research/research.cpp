/*
 * research.cpp
 *
 *  Created on: May 18, 2018
 *      Author: Naing Lin Tun
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <set>

using namespace std;

// function prototypes
vector<int> calcSquareMod(int); // return residue of square mods
int** adjMatrix(int, int, vector<int>, int&); // return pointer of adjacency matrix
void printAdjMatrix(int**, int num_nodes); // print results of adjacency matrix
vector<int> sievePrimes(int); // prime numbers up to int num
bool hamiltonian_cycle(int, int**, int, bool[], int);

void test_alg_prime(); // test conjecture for primes up to int num
void test_alg_nonprime(); // test conjecture from 2 to num
void test_alg_casebycase(); // test for each # of nodes and mod value

int main() {

	//test_alg_prime();
	test_alg_nonprime();
	//test_alg_casebycase();
}

vector<int> calcSquareMod(int mod_num) {
	/** This algorithm calculates the residue of square % mod_num excluding 0
	 *  For instance: if mod_num = 5:
	 *  	the residue of squares are mod of 1, 4, 9 16, 25, ....
	 *  	the residue set would be {1, 4} excluding 0
	 *
	 * 	@param: int mod_num
	 * 	@return: vector of residue set of squares
	 *
	 */

	bool residue[mod_num];
	vector<int> squareMod; // the residue set square mod

	for (int i = 1; i < mod_num; i++) {
		residue[i] = 0;
	}

	int x1=0, x2=0;
	int count = 0;
	for (int i = 1; i <= mod_num/2; i++) {
		int val = (i*i)%mod_num;
		x1 = val;

		if (x1 == x2 || (count >= 2 && val == squareMod[count-2])) {
			// The condition is for optimizing the algorithm.
			// Since the residue set (in modular arithmetic, Z/nZ) has a cyclic pattern,
			// we can break the loop if we detect "repetitions".
			// Repetition can occur in the following two forms:
			// 1. (expressed by the first disjunct of the condition):
			//		when two consecutive residue are equal
			//		for instance: (1^2) mod 3 = (2^2) mod 3.
			//		Then, the residue set is non other than {1} only.
			// 2. (expressed by the second disjunct of the condition):
			//		when the residue value "goes backward"
			// 		for instance, square mod 7 = {1, 4, 2, 4, 1, 4, 2, ....}
			//		We can see the cyclic pattern. Thus, the residue set is just
			//		{1, 4, 2}
			break;
		}

		if (residue[val] != 1 && val != 0) {
			// if the residue value is not already in the set (to avoid repetition)
			// and excluding 0 from the residue set
			residue[val] = 1;
			x2 = val;
			squareMod.push_back(val);
			count++;
		}
	}

	return squareMod;
}

int** adjMatrix(int num_nodes, int mod, vector<int> squareMod, int& num_edges) {
	/** This algorithm calculates the adjacency matrix for the following undirected graph:
	 * 		- there are 1, 2, ..., num_nodes vertices.
	 * 		- connect two vertices if we satisfy the following equation:
	 * 				sum of two vertices = square (% mod)
	 * 			or	(sum of two vertices % mod) = (square % mod)
	 * 		  In other words, if (sum of two vertices % mod) is one of the values
	 * 		  in of squareMod, then connect the two vertices via an edge.
	 *
	 *  @param: int num_nodes, int mod, vector of squareMod (the residue set), int& num_edges
	 *  	Note: num_edges is for tracking the number of edges. num_edges has to be
	 *  		defined before function call and pass it as a function parameter.
	 *  		Since it is passed by reference, the original num_edges will be changed
	 *  		as opposed to pass-by-value, which creates a copy of num_edges and the
	 *  		original num_edges remains unaffected.
	 * 	@return: pointer to the dynamically allocated 2d array,
	 * 			 which is the adjacency matrix of the graph
	 *
	 */

	// pointer to a dynamically allocated array, which consists of
	// pointers to the arrays too. It is an array of the array.
	// Thus, this is the 2d array

	int** adjMat = new int*[num_nodes];

	for (int i = 0; i < num_nodes; i++) { // initialize all entries of adjMat to 0
		adjMat[i] = new int[num_nodes];
		for (int j = 0; j < num_nodes; j++) {
			adjMat[i][j] = 0;
		}
	}

	for (int node_val = 1; node_val <= num_nodes; node_val++) {

		for (size_t i = 0; i < squareMod.size(); i++) {

			int connect_to = squareMod[i] - node_val;

			while (connect_to <= node_val) {
				// consider only edges where node_val > connect_to
				// to avoid repetition. For instance, 1--2 edge is
				// also the same as 2--1 edge; hence, redundancy.
				connect_to += mod;
			}

			while(connect_to <= num_nodes) {
				num_edges++;

				adjMat[node_val-1][connect_to-1] = 1; // -1 since index starts at 0
				adjMat[connect_to-1][node_val-1] = 1; // because adjMat is reflective diagonally

				connect_to += mod;
			}

		}

	}

	return adjMat;

}

void printAdjMatrix(int** ptr_adjMat, int num_nodes) {
	/** This algorithm prints the adjacency matrix, the neighborhood of each vertex
	 *  (both with and without repetitions), and the degree of each vertex.
	 *
	 *  @param: pointer to adjaceny matrix, the number of nodes
	 *  @return: print results as stated above
	 *
	 */


	cout << "Adjacency matrix: " << endl;

	for (int row = 0; row < num_nodes; row++) {
		for (int col = 0; col < num_nodes; col++) {
			cout << ptr_adjMat[row][col] << " ";
		}
		cout << endl;
	}

	cout << endl;

	for (int row = 0; row < num_nodes; row++) {

		int degree=0;
		cout << row+1 << " conects to ";
		for (int col = 0; col < num_nodes; col++) {
			if (ptr_adjMat[row][col] == 1) {
				degree++;
				cout << col+1 << " ";
			}
		}
		cout << "(Degree: " << degree << ")" << endl;

	}

}

vector<int> sievePrimes(int num_nodes) {
	/** This algorithm calculates the number of primes from 1 to num_nodes
	 *  using Sieve of Eratosthenes, with complexity O( n log(log n) ).
	 *
	 *  @param: int num_nodes
	 * 	@return: vector of primes from 1 to num_nodes, inclusively
	 */

	bool primeArray[num_nodes+1]; // declare array with type bool with size 1000

	// count the number of primes to store in vector<int> primes
	// -1 because 1 is not a prime
	int count = num_nodes-1;

	for (int i = 2; i <= num_nodes; ++i) {
		// initialize primeArray values from index 2 to 100 as true
		primeArray[i] = true;
	}


	for (size_t j = 2; j <= sqrt(num_nodes); ++j) {

		if (primeArray[j]) { // optimization: avoid repetition of visiting non-prime indexes
			for (size_t factorJ=2; factorJ <= num_nodes/j; factorJ++) {
				if (primeArray[factorJ * j]) {
					primeArray[factorJ * j] = false;
					count--;
				}
			}
		}

	}

	vector<int> primes(count);
	int index = 0;
	for (int i = 2; i <= num_nodes; ++i) {
		if (primeArray[i]) {
			primes[index] = i;
			index++;
		}
	}
	return primes;

}


bool hamiltonian_cycle(int start_node, int** ptr_adjMat, int num_nodes,
		bool visited[], int edge_count) {
	/** This algorithm checks whether a Hamiltonian cycle is present in a graph
	 *
	 *  @param: the start node, the adjacency matrix of the graph,
	 *  		the total number of nodes, an array to track
	 *  		visited vertices, edges counter
	 *
	 *  @return: true if Hamiltonian cycle is present; false otherwise
	 *
	 */

	if (edge_count == num_nodes-1 && ptr_adjMat[start_node][0] == 1 && num_nodes > 2) {
		// this is the base case.
		// 1. The first disjunct is true if we can find a Hamiltonian
		//		path from starting from start_node.
		// 2. The second disjunct is true if there is an edge
		//		connecting the start and end vertices of the path
		// 3. The third disjunct is for the case with just 1 or 2 vertices
		//		where Hamiltonian cycle cannot exist by definition
		//
		// Thus, the first conjunct must be true first in order for
		// the second conjunct to be true.
		//
		// So, "edge_count == num_nodes-1" is placed as the first conjunct
		// to take advantage of the "short-circuiting" property of &&.

		return true;
	}

	visited[start_node] = true;

	for (int col = 0; col < num_nodes; col++) {

		if (visited[col] == 0 && ptr_adjMat[start_node][col] == 1) {

			visited[col] = true;

			if (hamiltonian_cycle(col, ptr_adjMat, num_nodes, visited, edge_count+1)) {
				return true;
			}

			// if hamiltonian_cycle is false,
			// backtracks it
			visited[col] = false;

		}

	}

	return false;

}

void test_alg_prime() {
	/** Test algorithm for primes nodes from 1 to num
	 *
	 *  @return: print ( primes \in [1, num] : # of edges )
	 */

	int num;
	cout << "Enter num: ";
	cin >> num;

	vector<int> primes = sievePrimes(num);

	for (int i : primes) {
		int num_edges = 0;
		vector<int> squareMod = calcSquareMod(i);

		bool visited[i];
		for (int j = 0; j < i; j++) {
			visited[j] = 0;
		}
		//int visited[i];
		int** ptr_adjMat = adjMatrix(i, i, squareMod, num_edges);
		cout << i << " : " << num_edges << boolalpha << " "
						<< hamiltonian_cycle(0, ptr_adjMat, i, visited, 0)
						<< endl;

		delete[] ptr_adjMat;
		delete[] adjMatrix(i, i, squareMod, num_edges);
	}

}

void test_alg_nonprime() {
	/** Test algorithm for nodes from 1 to num
	 *
	 *  @return: print ( nodes \in [1, num] : # of edges )
	 */

	int num;
	cout << "Enter num: ";
	cin >> num;

	for (int i = 2; i <= num; i++) {
		int num_edges = 0;
		vector<int> squareMod = calcSquareMod(i);

		bool visited[i];
		for (int j = 0; j < i; j++) {
			visited[j] = 0;
		}
		//int visited[i];
		int** ptr_adjMat = adjMatrix(i, i, squareMod, num_edges);
		cout << i << " : " << num_edges << boolalpha << " "
				<< hamiltonian_cycle(0, ptr_adjMat, i, visited, 0)
				<< endl;

		//hamCycle(ptr_adjMat, i);

		delete[] ptr_adjMat;
		delete[] adjMatrix(i, i, squareMod, num_edges);
	}

}

void test_alg_casebycase() {
	/** Test algorithm for a specific case of num_nodes and mod
	 * 	For instance:
	 * 		number of nodes = 7, mod = 5.
	 * 		- Connect two nodes if (sum of two nodes % 5) is in
	 * 			the residue set of (square % 5)
	 *
	 *  @return: call printAdjMatrix() to print the graph's data
	 *  		and the number of edges
	 *
	 */

	int num_nodes, mod;
	cout << "Enter num_nodes and mod value, each separated by a space: ";
	cin >> num_nodes >> mod;
	cout << endl;

	vector<int> squareMod = calcSquareMod(mod);

	cout << "residue for square mod " << mod << ": ";
	for (int i : squareMod) {
		cout << i << " ";
	}

	cout << endl << endl;
	int num_edges = 0;

	int** ptr_adjMat = adjMatrix(num_nodes, mod, squareMod, num_edges);

	printAdjMatrix(ptr_adjMat, num_nodes);
	cout << "\nnum of edges: " << num_edges;


	bool visited[num_nodes];
	for (int i = 0; i < num_nodes; i++) {
		visited[i] = 0;
	}
	cout << endl;
	cout << boolalpha << "\nHamiltonian cycle present: "
			<< hamiltonian_cycle(0, ptr_adjMat, num_nodes, visited, 0) << endl << endl;

	delete[] ptr_adjMat;
	delete[] adjMatrix(num_nodes, mod, squareMod, num_edges);

}

