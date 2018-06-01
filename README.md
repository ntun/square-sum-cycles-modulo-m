# Computationally Finding Square-Sum Cycles modulo m

This program is written in C++ for the purpose of helping with the following problem:

Given a sequence of integers from {1, 2, 3, ..., n}, can we arrange the integers in a way that the (sum of two consecutive integers % m) equals (square % m), where square = {1, 4, 9, 16, 25, ...}? The first and last terms of the sequence also count as two consecutive terms. Mathematically, we are trying to solve the following equation:

<p align="center">
   <img src="http://latex.codecogs.com/gif.latex?x%20&plus;%20y%20%5Cequiv%20a%5E2%20%5Ctext%7B%20%28mod%20%7D%20m%20%5Ctext%7B%29%7D"> </p>

 where x, y are any two vertices in {1, 2, 3, ..., n} and a, m are positive integers.
 
 Our approach is to solve the problem using graph theory. For a given m and n, construct a graph, denoted ![equation](http://latex.codecogs.com/gif.latex?G_m%20%28n%29) with each term in a sequence as a vertex. Then, draw an edge between two vertices if and only if the above equation is satisfied. Next, we find the Hamiltonian cycle in the graph, and we are done.
 
 To give credit what's due, the motivation for this project is drawn from Numberphile [The Square-Sum Problem](https://www.youtube.com/watch?v=G1m7goLCJDY) video.
 
 # Getting Started
 
 This program is written in C++11 with MINGW compiler on a Windows-10 64-bit computer.
 
 Details for getting started with MINGW can be found [here](http://www.mingw.org/wiki/Getting_Started). You can download MINGW [here](https://sourceforge.net/projects/mingw/files/latest/download).
 
 I also use Eclipse program to write, test, debug, and develop C++ codes. For detailed information about getting started with Eclipse, see [here](https://www3.ntu.edu.sg/home/ehchua/programming/howto/EclipseCpp_HowTo.html) or [here](https://www.rose-hulman.edu/class/csse/resources/Eclipse/eclipse-c-configuration.htm).
 
 # How to use the program
 The program includes three test algorithms for three different cases:
 ```
 > test_alg_prime()
 > test_alg_nonprime()
 > test_alg_casebycase()
 ```
 
 Use the test algorithm only **one at a time**. Comment out the test algorithm in the *int main()* function for the ones that are not currently in use.
 ## Using test_alg_prime()
 
 This algorithm prompts the user for integer input first, call it num. Then, it checks whether a Hamiltonian cycle exists in ![equation](http://latex.codecogs.com/gif.latex?G_k%20%28k%29) for all k, where k is a prime between 1 and num.
 
 ## Using test_alg_nonprime()
 
 The same as test_alg_prime(), but k is now **all** integers between 2 and num inclusively.
 
 ## Using test_alg_casebycase()
 
 This algorithm prompts the user for number of nodes (n) and modulo value (m). Then, it checks whether a Hamiltonian cycle exists in ![equation](http://latex.codecogs.com/gif.latex?G_m%20%28n%29).
 
 # Authors
 - Naing Lin Tun
 
 # Acknowledgments
 - Professor Andrew Shallue from Illinois Wesleyan University for mentorship and advice on our research project in his MATH 370 course.
 - Patrick Ward, my project partner, for his contribution towards his mathematical findings and great ideas.
 - Isaac Dragomir for providing information about his previous talk on the similar topic "The Square-Sum Problem" 
